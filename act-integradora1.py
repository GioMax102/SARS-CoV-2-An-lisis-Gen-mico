def leer_secuencia(nombre):
    with open("data/" + nombre) as f:
        return "".join(l.strip() for l in f if not l.startswith(">")).upper()

def lps(P):
    m = len(P)
    tabla = [0] * m
    l, i = 0, 1
    while i < m:
        if P[i] == P[l]:
            l += 1
            tabla[i] = l
            i += 1
        elif l != 0:
            l = tabla[l - 1]
        else:
            i += 1
    return tabla


def kmp_search(text, pattern):
    n, m = len(text), len(pattern)
    lps_arr = lps(pattern)
    matches = []
    i = j = 0
    while i < n:
        if text[i] == pattern[j]:
            i += 1
            j += 1
        elif j > 0:
            j = lps_arr[j - 1]
        else:
            i += 1
        if j == m:
            matches.append(i - j)
            j = lps_arr[j - 1]
    return matches

wuhan = leer_secuencia("SARS-COV-2-MN908947.3.txt")

genes = {
    "M": leer_secuencia("gen-M.txt"),
    "S": leer_secuencia("gen-S.txt"),
    "ORF1AB": leer_secuencia("gen-ORF1AB.txt"),
}

for nombre, gen in genes.items():
    pos = kmp_search(wuhan, gen)
    print(f"Gen {nombre}")
    print(f"  Primeros 12: {gen[:12]}")
    for p in pos:
        print(f"  Inicio: {p}  Fin: {p + len(gen) - 1}")