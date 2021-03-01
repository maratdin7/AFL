#ifndef AFL_DATA_TYPES_H
#define AFL_DATA_TYPES_H

#include "glib.h"
#include "types.h"
#include "debug.h"

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


#define ptr_array_create(element_free_func) g_ptr_array_new_with_free_func(element_free_func)

#define ptr_array_insert(a, i, v) \
    g_ptr_array_insert(a, i, v)


#define ptr_array_set_size(a, s) g_ptr_array_set_size(a, s)

#define ptr_array_sort(a, compare) g_ptr_array_sort(a, compare)

#define ptr_array_remove_index(a, i) g_ptr_array_remove_index(a, i)

#define ptr_array_remove_index_fast(a, i) g_ptr_array_remove_index_fast(a, i)

#define ptr_array_index(a, i) g_ptr_array_index(a, i)

#define ptr_array_add(a, v) g_ptr_array_add(a, v)

#define ptr_array_set_index(a, i, val) do { \
    void *tmp = ptr_array_index(a, i); \
    *tmp = val; \
} while (0)

static u32 ptr_lower_bound(ptr_array_t *arr, u32 low, u32 hight, void *x, s32 (*compare)(void *, void *)) {

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

#define hash_table_get(ht, k) g_hash_table_lookup(ht, k)

static s32 int_equal(const void *a, const void *b) {
    return g_int_equal(a, b);
}

static u32 int_hash(const void *v) {
    return g_int_hash(v);
}

#define hash_table_create(hash_fun, equal_fun, key_free, val_free) \
        g_hash_table_new_full(hash_fun, equal_fun, key_free, val_free)

#define hash_table_lookup(hash_table, key) \
        g_hash_table_lookup(hash_table, key)

#define hash_table_lookup_extended(hash_table, key, orig_key, val) \
        g_hash_table_lookup_extended(hash_table, key, orig_key, val)

#define hash_table_destroy(ht) g_hash_table_destroy(ht)

#endif //AFL_DATA_TYPES_H
