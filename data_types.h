#ifndef AFL_DATA_TYPES_H
#define AFL_DATA_TYPES_H

#include "glib.h"

typedef GArray array_t;
typedef GPtrArray ptr_array_t;
typedef GHashTable hash_table_t;

#define array_index(a, t, i) ((a)->len > (i) ? &g_array_index(a, t, i) : NULL)


#define array_set_index(a, t, i, val) do { \
    t *_tmp = array_index(a, t, i); \
    *_tmp = val; \
} while (0)

#define array_set_size(a, s) g_array_set_size(a, s)

#define array_create(z, c, t) g_array_new(z,c, sizeof(t))

#define array_unref(a) g_array_unref(a)


#define ptr_array_create(fun) g_ptr_array_new_with_free_func(element_free_func)

#define ptr_array_insert(a, i, v) g_ptr_array_insert(a, i, v)

#define ptr_array_set_size(a, s) g_ptr_array_set_size(a, s)

void ptr_array_sort(ptr_array_t *a, s32 (compare) (const void *p1, const void *p2) ) {
    g_ptr_array_sort(a, compare);
}
array_t *ptr_array_remove_index(ptr_array_t *a, u32 i) {
    if (a->len <= i)
        FATAL("Out of bounds of array");
    return g_ptr_array_remove_index(a, i);
}

void *ptr_array_index(ptr_array_t *a, u32 i) {
    if (a->len <= i)
        FATAL("Out of bounds of array");

    return g_ptr_array_index(a, i);
}

#define ptr_array_set_index(a, i, val) do { \
    void *tmp = ptr_array_index(a, i); \
    *tmp = val; \
} while (0)

uint32_t ptr_lower_bound(ptr_array_t *arr, uint32_t low, uint32_t hight, void *x, int8_t (*compare)(void *, void *)) {

    while (low < hight) {

        uint32_t m = low + (hight - low) / 2;
        if (compare(x, ptr_array_index(arr, m)) <= 0)
            hight = m;
        else
            low = m + 1;
    }
    return low;
}

#define hash_table_insert(ht, k, v) g_hash_table_insert(ht, k, v)

#define hash_table_get(ht, k) g_hash_table_lookup(ht, k);

u32 double_hash(const void *d) {
    return g_double_hash(d);
}

s32 int_equal(const void *a, const void *b) {
    return g_int_equal(a, b);
}


#define hash_table_create(hash_fun, equal_fun, key_free, val_free) \
        g_hash_table_new_full(hash_fun, equal_fun, key_free, val_free)

#define hash_table_destroy(ht) g_hash_table_destroy(ht)

#endif //AFL_DATA_TYPES_H
