/* lmap_io_curl -- an HTTP(S) backend over libcurl range requests.  Optional: compiled
 * with -DLMAP_HTTP (make HTTP=1) and linked with -lcurl; otherwise lmap_io_open_url is a
 * stub that returns NULL, and the core library depends on nothing but zlib.
 *
 * Applications that already have a remote stack (htslib's hFILE, UCSC's udc) can plug
 * it in with lmap_io_open_backend instead; this one is for liftmap's own tools and for
 * callers that have none.
 */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#include "lmap_io.h"

#include <errno.h>
#include <stdlib.h>
#include <string.h>

#ifndef LMAP_HTTP

lmap_io *lmap_io_open_url(const char *url) { (void)url; errno = ENOSYS; return NULL; }
int lmap_http_available(void) { return 0; }

#else

#include <curl/curl.h>
#include <pthread.h>
#include <stdio.h>

typedef struct { CURL *h; pthread_mutex_t mu; } curl_ctx;
typedef struct { uint8_t *dst; int64_t want, got; } sink;

static size_t on_data(char *p, size_t size, size_t nmemb, void *ud) {
    sink *s = ud;
    size_t n = size * nmemb;
    if ((int64_t)n > s->want - s->got) return 0;        /* more than asked for: abort */
    memcpy(s->dst + s->got, p, n);
    s->got += (int64_t)n;
    return n;
}

static int64_t curl_read(void *vctx, void *buf, int64_t off, int64_t len) {
    curl_ctx *c = vctx;
    char range[64];
    snprintf(range, sizeof range, "%lld-%lld", (long long)off, (long long)(off + len - 1));
    sink s = { buf, len, 0 };
    pthread_mutex_lock(&c->mu);                          /* one easy handle, reused */
    curl_easy_setopt(c->h, CURLOPT_NOBODY, 0L);
    curl_easy_setopt(c->h, CURLOPT_HTTPGET, 1L);
    curl_easy_setopt(c->h, CURLOPT_RANGE, range);
    curl_easy_setopt(c->h, CURLOPT_WRITEFUNCTION, on_data);
    curl_easy_setopt(c->h, CURLOPT_WRITEDATA, &s);
    CURLcode rc = curl_easy_perform(c->h);
    long code = 0;
    curl_easy_getinfo(c->h, CURLINFO_RESPONSE_CODE, &code);
    pthread_mutex_unlock(&c->mu);
    /* 206 is a range answer; a server that ignores Range sends 200 and the whole body,
     * which on_data refuses unless the request was for all of it */
    if (rc != CURLE_OK || (code != 206 && !(code == 200 && off == 0))) return -1;
    return s.got;
}

static void curl_close(void *vctx) {
    curl_ctx *c = vctx;
    curl_easy_cleanup(c->h);
    pthread_mutex_destroy(&c->mu);
    free(c);
}

static const lmap_io_backend CURL_BE = { "http", curl_read, curl_close };

static pthread_once_t once = PTHREAD_ONCE_INIT;
static void global_init(void) { curl_global_init(CURL_GLOBAL_DEFAULT); }

lmap_io *lmap_io_open_url(const char *url) {
    pthread_once(&once, global_init);
    curl_ctx *c = calloc(1, sizeof *c);
    if (!c) return NULL;
    if (!(c->h = curl_easy_init()) || pthread_mutex_init(&c->mu, NULL) != 0) {
        if (c->h) curl_easy_cleanup(c->h);
        free(c); return NULL;
    }
    curl_easy_setopt(c->h, CURLOPT_URL, url);
    curl_easy_setopt(c->h, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(c->h, CURLOPT_FAILONERROR, 1L);
    curl_easy_setopt(c->h, CURLOPT_NOSIGNAL, 1L);        /* safe with threads */
    /* the size, from a HEAD request */
    curl_easy_setopt(c->h, CURLOPT_NOBODY, 1L);
    curl_off_t size = -1;
    if (curl_easy_perform(c->h) != CURLE_OK ||
        curl_easy_getinfo(c->h, CURLINFO_CONTENT_LENGTH_DOWNLOAD_T, &size) != CURLE_OK || size < 0) {
        curl_close(c); return NULL;
    }
    lmap_io *io = lmap_io_open_backend(&CURL_BE, c, (int64_t)size);
    if (!io) curl_close(c);
    return io;
}

int lmap_http_available(void) { return 1; }

#endif
