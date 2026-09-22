/* sections -- application sections and member-name prefix ranges. */
#include "../src/lmap_file.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int fails = 0;
#define CHECK(what, cond) do { if (!(cond)) { printf("  FAIL %s\n", what); fails++; } \
                               else printf("  ok   %s\n", what); } while (0)

int main(int argc, char **argv) {
    const char *path = argc > 1 ? argv[1] : "/tmp/lmap_sections_test.lmap";
    lmap_writer *w = lmap_writer_open(path, "test", "liftmap 1\n");
    const char *names[] = { "g2.x", "g1.chr2", "g10.chrA", "g1.chr1", "h.y" };
    for (int i = 0; i < 5; i++) lmap_writer_add_member(w, 0, names[i], 1000);
    lmap_writer_add_member(w, 1, "column", 100000);
    for (int i = 0; i < 5; i++) lmap_writer_add_run(w, i, 0, 10, 100 * i, 5, 0);

    const char meta[] = "T\t100000\nmax_gap\t0\n";
    uint8_t bin[300]; for (int i = 0; i < 300; i++) bin[i] = (uint8_t)(i * 7);
    printf("writer:\n");
    CHECK("adds a text section",           lmap_writer_add_section(w, "x.tui.meta", meta, strlen(meta)) == 0);
    CHECK("adds a binary section",         lmap_writer_add_section(w, "x.tui.anchor", bin, sizeof bin) == 0);
    CHECK("adds an empty section",         lmap_writer_add_section(w, "x.empty", NULL, 0) == 0);
    CHECK("refuses a duplicate id",        lmap_writer_add_section(w, "x.tui.meta", "z", 1) != 0);
    CHECK("refuses an id outside x.*",     lmap_writer_add_section(w, "dir.a", "z", 1) != 0);
    CHECK("closes",                        lmap_writer_close(w) == 0);

    lmap_file *f = lmap_open(path);
    printf("reader:\n");
    CHECK("opens", f != NULL);
    if (!f) return 1;
    uint8_t *p; size_t n;
    CHECK("reads the text section back",   lmap_read_section(f, "x.tui.meta", &p, &n) == 0 &&
                                           n == strlen(meta) && !memcmp(p, meta, n));
    if (n == strlen(meta)) free(p);
    CHECK("reads the binary section back", lmap_read_section(f, "x.tui.anchor", &p, &n) == 0 &&
                                           n == sizeof bin && !memcmp(p, bin, n));
    free(p);
    CHECK("reads the empty section",       lmap_read_section(f, "x.empty", &p, &n) == 0 && n == 0);
    free(p);
    CHECK("reports an absent section",     lmap_read_section(f, "x.nope", &p, &n) == 1);

    uint32_t r0, cnt;
    CHECK("prefix g1. = g1.chr1, g1.chr2 only (not g10.*)",
          lmap_member_prefix(f, 0, "g1.", &r0, &cnt) == 0 && cnt == 2 &&
          !strcmp(lmap_member_at(f, 0, lmap_member_by_rank(f, 0, r0))->name, "g1.chr1") &&
          !strcmp(lmap_member_at(f, 0, lmap_member_by_rank(f, 0, r0 + 1))->name, "g1.chr2"));
    CHECK("prefix g10. = g10.chrA",
          lmap_member_prefix(f, 0, "g10.", &r0, &cnt) == 0 && cnt == 1 &&
          !strcmp(lmap_member_at(f, 0, lmap_member_by_rank(f, 0, r0))->name, "g10.chrA"));
    CHECK("absent prefix gives an empty range", lmap_member_prefix(f, 0, "zz.", &r0, &cnt) == 0 && cnt == 0);
    CHECK("empty prefix covers every member",   lmap_member_prefix(f, 0, "", &r0, &cnt) == 0 && r0 == 0 && cnt == 5);
    CHECK("out-of-range rank refused",          lmap_member_by_rank(f, 0, 5) == -1);
    lmap_close(f);

    /* corrupt one byte inside the binary section: the section must fail its checksum */
    FILE *fp = fopen(path, "r+b");
    fseek(fp, 0, SEEK_END); long sz = ftell(fp); rewind(fp);
    uint8_t *all = malloc((size_t)sz);
    if (fread(all, 1, (size_t)sz, fp) != (size_t)sz) return 1;
    uint8_t *hit = memmem(all, (size_t)sz, bin, sizeof bin);
    if (hit) { hit[150] ^= 0xFF; rewind(fp); fwrite(all, 1, (size_t)sz, fp); }
    fclose(fp); free(all);
    f = lmap_open(path);
    printf("corruption:\n");
    CHECK("file with a corrupt app section still opens (loaded lazily)", f != NULL);
    CHECK("the corrupt section fails its checksum on read",
          f && hit && lmap_read_section(f, "x.tui.anchor", &p, &n) == -1);
    CHECK("an intact section still reads", f && lmap_read_section(f, "x.tui.meta", &p, &n) == 0);
    if (f) { free(p); lmap_close(f); }

    /* writer contract: runs must lie inside both members and not overlap on axis a */
    printf("run contract:\n");
    w = lmap_writer_open(path, "test", "liftmap 1\n");
    lmap_writer_add_member(w, 0, "s", 100);
    lmap_writer_add_member(w, 1, "column", 1000);
    CHECK("refuses a run past the axis-a member end", lmap_writer_add_run(w, 0, 0, 95, 0, 10, 0) != 0);
    CHECK("refuses a run past the axis-b member end", lmap_writer_add_run(w, 0, 0, 0, 995, 10, 1) != 0);
    CHECK("accepts a run ending exactly at both ends", lmap_writer_add_run(w, 0, 0, 90, 990, 10, 0) == 0);
    CHECK("refuses a run overlapping the last on axis a", lmap_writer_add_run(w, 0, 0, 95, 0, 1, 0) != 0);
    lmap_writer_abort(w);

    /* A_OVERLAP: axis a may overlap, order b only */
    printf("axis-a overlap:\n");
    w = lmap_writer_open(path, "test", "liftmap 1\n");
    CHECK("refused under order a",  lmap_writer_set_a_overlap(w) != 0);
    lmap_writer_set_order(w, LMAP_ORDER_B);
    CHECK("accepted under order b", lmap_writer_set_a_overlap(w) == 0);
    CHECK("order a is then refused", lmap_writer_set_order(w, LMAP_ORDER_A) != 0);
    lmap_writer_add_member(w, 0, "s", 100);
    lmap_writer_add_member(w, 1, "column", 1000);
    CHECK("an overlapping run is accepted",
          lmap_writer_add_run(w, 0, 0, 10, 500, 20, 0) == 0 &&
          lmap_writer_add_run(w, 0, 0, 15, 100, 10, 1) == 0 &&
          lmap_writer_add_run(w, 0, 0,  5, 800, 3, 0) == 0);
    CHECK("closes", lmap_writer_close(w) == 0);
    f = lmap_open(path);
    CHECK("reopens and reports the overlap", f && lmap_a_overlap(f) == 1);
    lmap_hit *hh = NULL; size_t hn = 0;
    CHECK("query_a returns every overlapping run, sorted by (a, b)",
          f && lmap_query_a(f, 0, 0, 100, &hh, &hn, NULL) == 0 && hn == 3 &&
          hh[0].a == 5 && hh[1].a == 10 && hh[1].b == 500 && hh[2].a == 15 && hh[2].strand == 1);
    free(hh);
    CHECK("query_a at a doubly-covered base returns both",
          f && lmap_query_a(f, 0, 17, 18, &hh, &hn, NULL) == 0 && hn == 2 &&
          hh[0].b == 107 && hh[1].b == 507);   /* tie on a: by b */
    free(hh);
    if (f) lmap_close(f);

    printf("%s (%d failures)\n", fails ? "FAILURES" : "ALL PASS", fails);
    return fails ? 1 : 0;
}
