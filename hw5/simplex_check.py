
import numpy as np
from itertools import product


def build_hamming_H(r):
    n = 2**r - 1
    H = np.zeros((r, n), dtype=int)
    col = 0
    for i in range(1, 2**r):
        bits = []
        val = i
        for _ in range(r):
            bits.append(val % 2)
            val //= 2
        H[:, col] = bits[::-1]
        col += 1
    return H


def generate_all_codewords(G):
    k, n = G.shape
    codewords = []
    for m in product([0, 1], repeat=k):
        m_vec = np.array(m, dtype=int)
        c = m_vec @ G % 2
        codewords.append((m_vec, c))
    return codewords


def hamming_weight(c):
    return int(np.sum(c))


def hamming_distance(c1, c2):
    return int(np.sum((c1 + c2) % 2))


def check_simplex(r):
    print("=" * 70)
    print(f"  Параметр r = {r}")
    print(f"  Код Хэмминга: ({2**r - 1}, {2**r - 1 - r}), d_min = 3")
    print(f"  Дуальный код: ({2**r - 1}, {r})")
    print(f"  Ожидаемый вес ненулевых слов: 2^(r-1) = {2**(r-1)}")
    print("=" * 70)

    H = build_hamming_H(r)
    n = 2**r - 1
    k = r

    print(f"\nПроверочная матрица H кода Хэмминга ({r} x {n}):")
    for i in range(r):
        print(f"  h{i+1}: {[int(j) for j in H[i]]}")

    G_dual = H.copy()
    print(f"\nПорождающая матрица дуального кода G_dual = H ({k} x {n}):")
    for i in range(k):
        print(f"  g{i+1}: {[int(j) for j in G_dual[i]]}")

    codewords = generate_all_codewords(G_dual)
    M = len(codewords)
    print(f"\nВсего кодовых слов: {M} = 2^{k}")

    print(f"\n{'m':>12s}  |  {'кодовое слово':>25s}  |  вес")
    print("-" * 55)

    weights = []
    nonzero_weights = []
    for m_vec, c in codewords:
        w = hamming_weight(c)
        weights.append(w)
        m_str = ''.join(map(str, m_vec))
        c_str = ''.join(map(str, c))
        print(f"  {m_str:>10s}  |  {c_str:>25s}  |  {w}")
        if w > 0:
            nonzero_weights.append(w)

    expected_weight = 2**(r - 1)
    all_same = all(w == expected_weight for w in nonzero_weights)

    print(f"\n--- Результаты проверки ---")
    print(f"  Ожидаемый вес ненулевых слов: {expected_weight}")
    print(f"  Фактические веса ненулевых слов: {sorted(set(nonzero_weights))}")
    print(f"  Все ненулевые слова имеют вес {expected_weight}: "
          f"{'ДА' if all_same else 'НЕТ'}")

    print(f"\n--- Проверка расстояний ---")
    distances = set()
    for i in range(M):
        for j in range(i + 1, M):
            _, ci = codewords[i]
            _, cj = codewords[j]
            d = hamming_distance(ci, cj)
            if d > 0:
                distances.add(d)

    print(f"  Множество расстояний между различными словами: {sorted(distances)}")

    all_equal_dist = (len(distances) == 1)
    print(f"  Все расстояния одинаковы: {'ДА' if all_equal_dist else 'НЕТ'}")

    print(f"\n--- Распределение весов ---")
    from collections import Counter
    weight_counts = Counter(weights)
    for w in sorted(weight_counts.keys()):
        print(f"  Вес {w}: {weight_counts[w]} слов")

    is_simplex = all_same and all_equal_dist
    print(f"\n{'=' * 70}")
    if is_simplex:
        print(f"  ВЫВОД: Дуальный код к коду Хэмминга (r={r}) ЯВЛЯЕТСЯ симплексным.")
        print(f"  Все {len(nonzero_weights)} ненулевых слов имеют вес {expected_weight}.")
        print(f"  Расстояние между любыми двумя различными словами = {expected_weight}.")
    else:
        print(f"  ВЫВОД: Код НЕ является симплексным.")
    print(f"{'=' * 70}\n")

    return is_simplex


if __name__ == "__main__":
    print("ПРОВЕРКА: дуальный код к коду Хэмминга — симплексный код\n")

    results = {}
    for r in [2, 3, 4, 5, 6, 7]:
        results[r] = check_simplex(r)

    print("\n" + "=" * 70)
    print("  ИТОГОВАЯ ТАБЛИЦА")
    print("=" * 70)
    print(f"  {'r':>3s}  {'(n, k)':>10s}  {'d_min':>6s}  {'Вес':>6s}  {'Симплекс?':>10s}")
    print("-" * 50)
    for r in [2, 3, 4, 5]:
        n = 2**r - 1
        k = r
        w = 2**(r - 1)
        ok = "ДА" if results[r] else "НЕТ"
        print(f"  {r:>3d}  ({n:>2d}, {k:>2d})  {w:>6d}  {w:>6d}  {ok:>10s}")
    print()
