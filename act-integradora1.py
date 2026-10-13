import os

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, "data")


# ===============================================================
# LECTURA DE ARCHIVOS
# Abre los archivos de la carpeta data/ y los convierte en cadenas
# o diccionarios para poder trabajar con ellos.
# ===============================================================

def leer_secuencia(nombre):
    with open(os.path.join(DATA, nombre), encoding="utf-8") as f:
        return "".join(
            l.strip() for l in f if not l.startswith(">")
        ).upper()


def leer_proteinas(nombre):
    proteinas = {}
    nombre_actual = None

    with open(os.path.join(DATA, nombre), encoding="utf-8") as archivo:
        for linea in archivo:
            linea = linea.strip()

            if not linea:
                continue

            if linea.startswith(">"):
                nombre_actual = linea[1:].split()[0]
                proteinas[nombre_actual] = ""
            elif nombre_actual is not None:
                proteinas[nombre_actual] += linea.upper().replace(" ", "")

    return proteinas


# ===============================================================
# PUNTO 1 - BUSQUEDA DE GENES
# Busca cada gen (M, S y ORF1AB) dentro del genoma de Wuhan usando
# KMP, O(n + m). Muestra el nombre, los indices de aparicion y los
# primeros 12 caracteres.
# ===============================================================

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
    if m == 0 or m > n:
        return []

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


def mostrar_punto1(wuhan, genes):
    print("\n--- PUNTO 1: BUSQUEDA DE GENES ---")

    for nombre, gen in genes.items():
        posiciones = kmp_search(wuhan, gen)

        print(f"\nGen {nombre}")
        print(f"  Primeros 12: {gen[:12]}")

        if not posiciones:
            print("  No encontrado")
        for p in posiciones:
            print(f"  Inicio: {p}  Fin: {p + len(gen) - 1}")

# ===============================================================
# PUNTO 2 - PALINDROMO GENÉTICO MAS LARGO
# Encuentra el palíndromo genético más largo de cada gen con Manacher.
# Complejidad: O(n) tiempo, donde n es la longitud de la cadena del gen.
# ===============================================================

def manacher_genetico(cadena):
    complemento = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    texto = "@$" + "$".join(cadena) + "$#"
    e = len(texto)

    P = [0] * e
    centro = 0
    limite = 0

    def coinciden(izq, der):
        if izq == '$' and der == '$':
            return True
        if izq in complemento and der in complemento:
            return complemento[izq] == der
        return False

    for i in range(1, e - 1):
        if texto[i] != '$':
            continue

        if i < limite:
            simetrica = 2 * centro - i
            P[i] = min(limite - i, P[simetrica])

        gap = P[i] + 1

        while coinciden(texto[i - gap], texto[i + gap]):
            P[i] += 1
            gap += 1

        if i + P[i] > limite:
            limite = i + P[i]
            centro = i

    longitud_maxima = max(P)
    posicion = P.index(longitud_maxima)

    inicio = (posicion - longitud_maxima) // 2
    palindromo = cadena[inicio:inicio + longitud_maxima]

    return palindromo, longitud_maxima, inicio


def mostrar_punto2(genes):
    print("\n--- PUNTO 2: PALINDROMOS GENÉTICOS ---")
    BASE = os.path.dirname(os.path.abspath(__file__)) # Asegurar ruta base

    with open(os.path.join(BASE, "palindromos.txt"), "w", encoding="utf-8") as archivo:
        for nombre, gen in genes.items():
            palindromo, longitud, ini = manacher_genetico(gen)

            print(f"\nGen {nombre}")
            print(f"  Longitud: {longitud}")
            print(f"  Inicio dentro del gen: {ini}")
            print(f"  Palindromo: {palindromo}")

            archivo.write(f"Gen {nombre}\n")
            archivo.write(f"Longitud: {longitud}\n")
            archivo.write(f"Inicio dentro del gen: {ini}\n")
            archivo.write(f"Palindromo: {palindromo}\n\n")

    print("\nPalindromos guardados en palindromos.txt")


# ===============================================================
# PUNTO 3 - IDENTIFICACION DE PROTEINAS
# Traduce el genoma en los 3 marcos de lectura y busca cada proteina
# con KMP para saber de que parte del genoma se produce. La proteina
# nsp12 (QHD43415_11) tiene un corrimiento ribosomal -1, por eso se
# busca aparte. Muestra indices, primeros 4 aminoacidos y codones.
# ===============================================================

def crear_tabla_codones():
    bases = "TCAG"
    aminoacidos = (
        "FFLLSSSSYY**CC*W"
        "LLLLPPPPHHQQRRRR"
        "IIIMTTTTNNKKSSRR"
        "VVVVAAAADDEEGGGG"
    )

    tabla = {}
    indice = 0

    for primera in bases:
        for segunda in bases:
            for tercera in bases:
                tabla[primera + segunda + tercera] = aminoacidos[indice]
                indice += 1

    return tabla


def traducir_adn(adn, marco, tabla):
    proteina = []

    for i in range(marco, len(adn) - 2, 3):
        proteina.append(tabla.get(adn[i:i + 3], "X"))

    return "".join(proteina)


def buscar_con_corrimiento(genoma, secuencia, traducciones, tabla):
    largo = len(secuencia)
    if largo < 60:
        return None

    mitad = largo // 2
    ancla = secuencia[mitad:mitad + 30]

    for marco in range(3):
        for pos in kmp_search(traducciones[marco], ancla):
            k = mitad
            p = pos
            while k > 0 and p > 0 and \
                    traducciones[marco][p - 1] == secuencia[k - 1]:
                k -= 1
                p -= 1

            if k == 0:
                continue

            inicio_sufijo = marco + p * 3
            inicio_prefijo = inicio_sufijo + 1 - 3 * k
            if inicio_prefijo < 0:
                continue

            prefijo = "".join(
                tabla.get(genoma[i:i + 3], "X")
                for i in range(inicio_prefijo, inicio_prefijo + 3 * k, 3)
            )
            sufijo_ok = traducciones[marco][p:p + (largo - k)] == secuencia[k:]

            if prefijo == secuencia[:k] and sufijo_ok:
                fin = inicio_sufijo + (largo - k) * 3 - 1
                n_aa = min(4, largo)
                return {
                    "inicio": inicio_prefijo,
                    "fin": fin,
                    "marco": inicio_prefijo % 3,
                    "codones": [genoma[i:i + 3]
                                for i in range(inicio_prefijo,
                                               inicio_prefijo + n_aa * 3, 3)],
                    "corrimiento": (inicio_sufijo, k),
                }
    return None


def buscar_proteinas(genoma, proteinas, tabla):
    encontradas = {}
    no_encontradas = []

    traducciones = [traducir_adn(genoma, marco, tabla) for marco in range(3)]

    for nombre, secuencia in proteinas.items():
        coincidencias = []

        for marco in range(3):
            for posicion in kmp_search(traducciones[marco], secuencia):
                inicio = marco + posicion * 3
                fin = inicio + len(secuencia) * 3 - 1

                n_aa = min(4, len(secuencia))
                codones = [genoma[i:i + 3]
                           for i in range(inicio, inicio + n_aa * 3, 3)]

                coincidencias.append({
                    "inicio": inicio,
                    "fin": fin,
                    "marco": marco,
                    "codones": codones,
                })

        if not coincidencias:
            extra = buscar_con_corrimiento(
                genoma, secuencia, traducciones, tabla
            )
            if extra:
                coincidencias.append(extra)

        if coincidencias:
            encontradas[nombre] = coincidencias
        else:
            no_encontradas.append(nombre)

    return encontradas, no_encontradas


def mostrar_punto3(proteinas, encontradas, no_encontradas):
    print("\n--- PUNTO 3: PROTEINAS ---")

    for nombre, coincidencias in encontradas.items():
        print(f"\nProteina {nombre}")
        print(f"  Primeros 4 aminoacidos: {proteinas[nombre][:4]}")

        for r in coincidencias:
            print(f"  Inicio: {r['inicio']}")
            print(f"  Fin: {r['fin']}")
            print(f"  Marco: {r['marco']}")
            print(f"  Codones: {' '.join(r['codones'])}")
            if "corrimiento" in r:
                print(f"  Nota: corrimiento ribosomal -1 en el indice "
                      f"{r['corrimiento'][0]} (tras {r['corrimiento'][1]} aa)")

    print("\n--- RESUMEN PROTEINAS ---")
    print(f"Total de proteinas: {len(proteinas)}")
    print(f"Encontradas: {len(encontradas)}")
    print(f"No encontradas: {len(no_encontradas)}")

    if no_encontradas:
        print("\nProteinas no encontradas:")
        for nombre in no_encontradas:
            print(f"  {nombre}")


# ===============================================================
# PUNTO 4 - COMPARACION DE GENOMAS
# Compara Wuhan 2019 vs Texas 2020 posicion por posicion, O(n).
# Lista las diferencias, el tramo mas largo sin diferencias y, para
# cada diferencia, el codon afectado y el aminoacido de cada version
# (usando el marco de lectura de la proteina donde cae).
# ===============================================================

def comparar_genomas(wuhan, texas):
    diferencias = []
    inicio_actual = 0
    mejor_inicio = 0
    mejor_longitud = 0
    n = min(len(wuhan), len(texas))

    for i in range(n):
        if wuhan[i] != texas[i]:
            diferencias.append((i, wuhan[i], texas[i]))

            longitud = i - inicio_actual
            if longitud > mejor_longitud:
                mejor_longitud = longitud
                mejor_inicio = inicio_actual

            inicio_actual = i + 1

    longitud = n - inicio_actual
    if longitud > mejor_longitud:
        mejor_longitud = longitud
        mejor_inicio = inicio_actual

    return diferencias, mejor_inicio, mejor_longitud


def construir_regiones(encontradas):
    regiones = []
    for nombre, coincidencias in encontradas.items():
        for c in coincidencias:
            if "corrimiento" in c:
                inicio_sufijo, k = c["corrimiento"]
                regiones.append((c["inicio"], inicio_sufijo, nombre))
                regiones.append((inicio_sufijo, c["fin"], nombre))
            else:
                regiones.append((c["inicio"], c["fin"], nombre))
    return regiones


def buscar_region(i, regiones):
    mejor = None
    for inicio, fin, nombre in regiones:
        if inicio <= i <= fin:
            if mejor is None or (fin - inicio) < (mejor[1] - mejor[0]):
                mejor = (inicio, fin, nombre)
    return mejor


def analizar_mutaciones(wuhan, texas, diferencias, regiones, tabla):
    resultados = []

    for i, base_w, base_t in diferencias:
        region = buscar_region(i, regiones)

        if region is None:
            resultados.append({"indice": i, "base_w": base_w,
                               "base_t": base_t, "region": None})
            continue

        inicio, fin, nombre = region
        ini_codon = inicio + ((i - inicio) // 3) * 3

        cod_w = wuhan[ini_codon:ini_codon + 3]
        cod_t = texas[ini_codon:ini_codon + 3]
        aa_w = tabla.get(cod_w, "X")
        aa_t = tabla.get(cod_t, "X")

        resultados.append({
            "indice": i, "base_w": base_w, "base_t": base_t,
            "region": (inicio, fin, nombre), "ini_codon": ini_codon,
            "cod_w": cod_w, "cod_t": cod_t, "aa_w": aa_w, "aa_t": aa_t,
            "sinonima": aa_w == aa_t,
        })

    return resultados


def mostrar_punto4(wuhan, texas, diferencias, inicio, longitud):
    print("\n--- PUNTO 4: COMPARACION DE GENOMAS ---")

    print(f"\nLongitud Wuhan: {len(wuhan)}")
    print(f"Longitud Texas: {len(texas)}")

    if len(wuhan) != len(texas):
        print(f"NOTA: Texas es {len(wuhan) - len(texas)} nt mas corto; la "
              "diferencia esta en el extremo 3' (cola poly-A), por lo que se "
              "compara la region comun posicion a posicion.")

    print(f"\nDiferencias por posicion: {len(diferencias)}")

    with open(os.path.join(BASE, "diferencias.txt"), "w",
              encoding="utf-8") as archivo:
        for indice, base_w, base_t in diferencias:
            print(f"Indice {indice}: {base_w} -> {base_t}")
            archivo.write(f"{indice}\t{base_w}\t{base_t}\n")

    print("\nSegmento identico mas largo")
    print(f"Inicio: {inicio}")
    print(f"Fin: {inicio + longitud - 1}")
    print(f"Longitud: {longitud}")


def mostrar_mutaciones(resultados):
    print("\n--- ANALISIS DE MUTACIONES ---")

    sinonimas = 0
    no_sinonimas = 0
    fuera = 0

    for r in resultados:
        print(f"\nIndice: {r['indice']}  ({r['base_w']} -> {r['base_t']})")

        if r["region"] is None:
            fuera += 1
            print("  Fuera de las proteinas conocidas (region no codificante)")
            continue

        inicio, fin, nombre = r["region"]
        print(f"  Proteina: {nombre}  (region {inicio}-{fin})")
        print(f"  Posicion del codon: {r['ini_codon']}")
        print(f"  Wuhan: {r['cod_w']} -> {r['aa_w']}")
        print(f"  Texas: {r['cod_t']} -> {r['aa_t']}")

        if r["sinonima"]:
            sinonimas += 1
            print("  No cambia el aminoacido (mutacion sinonima)")
        else:
            no_sinonimas += 1
            print("  Cambia el aminoacido (mutacion no sinonima)")

    print("\n--- RESUMEN DE MUTACIONES ---")
    print(f"Sinonimas: {sinonimas}")
    print(f"No sinonimas: {no_sinonimas}")
    print(f"Fuera de proteinas conocidas: {fuera}")


# ===============================================================
# MAIN
# Lee los archivos y corre los 4 puntos en orden.
# ===============================================================

def main():
    wuhan = leer_secuencia("SARS-COV-2-MN908947.3.txt")
    texas = leer_secuencia("SARS-COV-2-MT106054.1.txt")

    genes = {
        "M": leer_secuencia("gen-M.txt"),
        "S": leer_secuencia("gen-S.txt"),
        "ORF1AB": leer_secuencia("gen-ORF1AB.txt"),
    }

    proteinas = leer_proteinas("seq-proteins.txt")
    tabla_codones = crear_tabla_codones()

    mostrar_punto1(wuhan, genes)

    mostrar_punto2(genes)

    encontradas, no_encontradas = buscar_proteinas(
        wuhan, proteinas, tabla_codones
    )
    mostrar_punto3(proteinas, encontradas, no_encontradas)

    diferencias, inicio, longitud = comparar_genomas(wuhan, texas)
    mostrar_punto4(wuhan, texas, diferencias, inicio, longitud)

    regiones = construir_regiones(encontradas)
    resultados = analizar_mutaciones(
        wuhan, texas, diferencias, regiones, tabla_codones
    )
    mostrar_mutaciones(resultados)


if __name__ == "__main__":
    main()