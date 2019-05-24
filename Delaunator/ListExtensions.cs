using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace Delaunator {
    internal static class ListExtensions {
        public static List<T> Fill<T>(this List<T> list, T value = default) {
            for (int i = 0; i < list.Capacity; i++) {
                list.Add(value);
            }
            return list;
        }
        public static List<T> Grow<T>(this List<T> list, int size) {
            if (size > list.Count) {
                if (list.Capacity < size) {
                    list.Capacity = size + (size / 2);
                }
                int count = size - list.Count;
                for (int i = 0; i < count; i++) {
                    list.Add(default(T));
                }
            }
            return list;
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void SetSafely<T>(this List<T> list, int index, T value) {
            if (index >= list.Count) {
                list.Grow(index + 1);
            }
            list[index] = value;
        }
    }
}