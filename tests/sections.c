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
    lmap_writer *w = lmap_writer_open(path, "test");
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
    w = lmap_writer_open(path, "test");
    lmap_writer_add_member(w, 0, "s", 100);
    lmap_writer_add_member(w, 1, "column", 1000);
    CHECK("refuses a run past the axis-a member end", lmap_writer_add_run(w, 0, 0, 95, 0, 10, 0) != 0);
    CHECK("refuses a run past the axis-b member end", lmap_writer_add_run(w, 0, 0, 0, 995, 10, 1) != 0);
    CHECK("accepts a run ending exactly at both ends", lmap_writer_add_run(w, 0, 0, 90, 990, 10, 0) == 0);
    CHECK("refuses a run overlapping the last on axis a", lmap_writer_add_run(w, 0, 0, 95, 0, 1, 0) != 0);
    lmap_writer_abort(w);

    /* A_OVERLAP: axis a may overlap, order b only */
    printf("axis-a overlap:\n");
    w = lmap_writer_open(path, "test");
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

    printf("names:\n");
    remove(path);
    w = lmap_writer_open(path, "test");
    lmap_writer_add_member(w, 0, "dup", 10); lmap_writer_add_member(w, 0, "dup", 10);
    lmap_writer_add_member(w, 1, "column", 100);
    CHECK("refuses duplicate member names at close", lmap_writer_close(w) != 0);
    { FILE *t = fopen(path, "rb"); char pp[4096]; snprintf(pp, sizeof pp, "%s.partial", path);
      FILE *q = fopen(pp, "rb");
      CHECK("a failed close leaves neither the file nor a partial", !t && !q);
      if (t) fclose(t);
      if (q) fclose(q); }

    printf("builder:\n");
    {
        lmap_builder *bd = lmap_builder_open(path, "test", 0, NULL);
        lmap_builder_add_member(bd, 0, "s", 100);
        lmap_builder_add_member(bd, 1, "t", 100);
        CHECK("refuses a run on an undeclared member", lmap_builder_add_run(bd, 0, 5, 0, 0, 10, 0) != 0 &&
              strstr(lmap_builder_error(bd), "not declared") != NULL);
        lmap_builder_abort(bd);
        bd = lmap_builder_open(path, "test", 0, NULL);
        lmap_builder_add_member(bd, 0, "s", 100);
        lmap_builder_add_member(bd, 1, "t", 100);
        CHECK("refuses a run past a member end", lmap_builder_add_run(bd, 0, 0, 95, 0, 10, 0) != 0 &&
              strstr(lmap_builder_error(bd), "outside") != NULL);
        lmap_builder_abort(bd);
        bd = lmap_builder_open(path, "test", 0, NULL);
        lmap_builder_add_member(bd, 0, "s", 100);
        lmap_builder_add_member(bd, 1, "t", 100);
        /* one stretch in pieces, out of order, with a duplicate: one canonical run */
        lmap_builder_add_run(bd, 0, 0, 20, 30, 5, 0);
        lmap_builder_add_run(bd, 0, 0, 10, 20, 10, 0);
        lmap_builder_add_run(bd, 0, 0, 12, 22, 3, 0);
        lmap_builder_add_run(bd, 0, 0, 50, 40, 5, 1);     /* reverse: 50..55 <-> 44..40 */
        lmap_builder_add_run(bd, 0, 0, 55, 35, 5, 1);     /* continues it: 55..60 <-> 39..35 */
        lmap_build_stats bs;
        CHECK("closes", lmap_builder_close(bd, &bs) == 0 && bs.error[0] == 0);
        CHECK("5 pieces -> 2 maximal runs, order a", bs.input_runs == 5 && bs.runs == 2 && bs.order == LMAP_ORDER_A);
        f = lmap_open(path);
        lmap_hit *hh; size_t hn;
        CHECK("the runs are [10,25)->[20,35) and [50,60)->[35,45) reversed",
              f && lmap_query_a(f, 0, 0, 100, &hh, &hn, NULL) == 0 && hn == 2 &&
              hh[0].a == 10 && hh[0].b == 20 && hh[0].len == 15 && hh[0].strand == 0 &&
              hh[1].a == 50 && hh[1].b == 35 && hh[1].len == 10 && hh[1].strand == 1);
        free(hh); lmap_close(f);
        bd = lmap_builder_open(path, "test", 0, NULL);
        lmap_builder_add_member(bd, 0, "s", 100);
        lmap_builder_add_member(bd, 1, "t", 100);
        lmap_builder_add_run(bd, 0, 0, 10, 20, 10, 0);
        lmap_builder_add_run(bd, 0, 0, 15, 60, 10, 0);
        CHECK("axis-a overlap fails close unless allowed, with a reason",
              lmap_builder_close(bd, &bs) != 0 && bs.overlapping_runs == 1 && strstr(bs.error, "overlap"));
        bd = lmap_builder_open(path, "test", 0, NULL);
        lmap_builder_add_member(bd, 0, "s", 100);
        lmap_builder_add_member(bd, 1, "t", 100);
        lmap_builder_allow_overlap(bd, 1);
        lmap_builder_add_run(bd, 0, 0, 10, 20, 10, 0);
        lmap_builder_add_run(bd, 0, 0, 15, 60, 10, 0);
        CHECK("allowed: written in order b with the overlap flag",
              lmap_builder_close(bd, &bs) == 0 && bs.order == LMAP_ORDER_B &&
              (f = lmap_open(path)) != NULL && lmap_a_overlap(f) == 1);
        if (f) lmap_close(f);
    }

    /* metadata and groups */
    printf("metadata and groups:\n");
    w = lmap_writer_open(path, "test");
    CHECK("sets a key",                     lmap_writer_set_meta(w, "tui.format", "0.4") == 0);
    CHECK("replaces a key",                 lmap_writer_set_meta(w, "tui.format", "0.5") == 0);
    CHECK("keeps tabs and newlines in values", lmap_writer_set_meta(w, "axis.a", "seq\tuence\nx") == 0);
    CHECK("refuses an empty key",           lmap_writer_set_meta(w, "", "v") != 0);
    lmap_writer_set_meta(w, "axis.b", "column");
    const char *seqs[] = { "gB.chr2", "gA.chr1", "gB.chr1", "loose", "gA.chr10" };
    const uint64_t lens[] = { 200, 100, 300, 50, 400 };
    for (int i = 0; i < 5; i++) lmap_writer_add_member(w, 0, seqs[i], lens[i]);
    lmap_writer_add_member(w, 1, "column", 100000);
    int32_t gB = lmap_writer_add_group(w, 0, "gB"), gA = lmap_writer_add_group(w, 0, "gA");
    CHECK("groups get ids in declaration order", gB == 0 && gA == 1);
    CHECK("refuses a duplicate group name", lmap_writer_add_group(w, 0, "gA") < 0);
    CHECK("refuses a member id out of range", lmap_writer_set_member_group(w, 0, 9, 0) != 0);
    CHECK("refuses a group id out of range",  lmap_writer_set_member_group(w, 0, 0, 7) != 0);
    lmap_writer_set_member_group(w, 0, 0, (uint32_t)gB); lmap_writer_set_member_group(w, 0, 2, (uint32_t)gB);
    lmap_writer_set_member_group(w, 0, 1, (uint32_t)gA); lmap_writer_set_member_group(w, 0, 4, (uint32_t)gA);
    for (int i = 0; i < 5; i++) lmap_writer_add_run(w, (uint32_t)i, 0, 0, 1000 * i, 10, 0);
    CHECK("closes", lmap_writer_close(w) == 0);
    f = lmap_open(path);
    CHECK("reopens", f != NULL);
    if (!f) return 1;
    CHECK("reads a key back (latest value)", lmap_meta_get(f, "tui.format") && !strcmp(lmap_meta_get(f, "tui.format"), "0.5"));
    CHECK("value keeps tab and newline", lmap_meta_get(f, "axis.a") && !strcmp(lmap_meta_get(f, "axis.a"), "seq\tuence\nx"));
    CHECK("absent key is NULL", lmap_meta_get(f, "nope") == NULL);
    const char *k0, *k2;
    CHECK("enumerates 3 keys in key order", lmap_meta_count(f) == 3 && lmap_meta_at(f, 0, &k0, NULL) == 0 &&
          lmap_meta_at(f, 2, &k2, NULL) == 0 && !strcmp(k0, "axis.a") && !strcmp(k2, "tui.format"));
    CHECK("2 groups on axis a, none on b", lmap_n_groups(f, 0) == 2 && lmap_n_groups(f, 1) == 0);
    CHECK("group names and lookup", !strcmp(lmap_group_name(f, 0, 0), "gB") && lmap_group_by_name(f, 0, "gA") == 1 &&
          lmap_group_by_name(f, 0, "gC") == -1);
    const uint32_t *mm; uint32_t nmm; uint64_t tl;
    CHECK("members in name order, with total length",
          lmap_group_members(f, 0, 1, &mm, &nmm, &tl) == 0 && nmm == 2 && tl == 500 &&
          !strcmp(lmap_member_at(f, 0, mm[0])->name, "gA.chr1") && !strcmp(lmap_member_at(f, 0, mm[1])->name, "gA.chr10"));
    CHECK("member -> group, and -1 for none",
          lmap_member_group(f, 0, (uint32_t)lmap_member_by_name(f, 0, "gB.chr1")) == 0 &&
          lmap_member_group(f, 0, (uint32_t)lmap_member_by_name(f, 0, "loose")) == -1);
    lmap_cursor *gc = lmap_cursor_open(f, 1, 0, mm, nmm, 0);
    lmap_hit gh; int ngh = 0; int only_a = 1;
    lmap_cursor_seek(gc, 0, 100000);
    while (lmap_cursor_next(gc, &gh) == 1) { ngh++; if (lmap_member_group(f, 0, gh.a_member) != 1) only_a = 0; }
    CHECK("a cursor restricted to a group sees only its members", ngh == 2 && only_a);
    lmap_cursor_close(gc);
    lmap_close(f);

    printf("%s (%d failures)\n", fails ? "FAILURES" : "ALL PASS", fails);
    return fails ? 1 : 0;
}
