# SARS-CoV-2 Genome Analysis

Análisis computacional del genoma del SARS-CoV-2 mediante algoritmos de procesamiento de strings.

## Descripción
Programa que analiza las secuencias genómicas de Wuhan (2019) y Texas (2020), localiza genes, detecta palíndromos genéticos, mapea proteínas y compara diferencias entre ambas versiones del virus.

## Tareas principales
- Localizar índices de los genes **M**, **S** y **ORF1ab** en el genoma.
- Encontrar el **palíndromo genético** (de complemento inverso) más largo en cada gen.
- Identificar secciones del genoma que producen cada proteína y sus codones.
- Comparar genomas de Wuhan 2019 vs. Texas 2020: diferencias, codones y aminoácidos afectados.

## Algoritmos
| Punto | Algoritmo | Complejidad |
|-------|-----------|-------------|
| 1. Búsqueda de genes | KMP | O(n + m) |
| 2. Palíndromo genético más largo | Manacher adaptado a complemento inverso | O(n) |
| 3. Identificación de proteínas | Traducción en 3 marcos de lectura + KMP | O(n + m) por marco |
| 4. Comparación de genomas | Comparación posición por posición | O(n) |

Todos los algoritmos están implementados desde cero, sin librerías que los resuelvan.

### Sobre los palíndromos (punto 2)
En biología molecular, un palíndromo de ADN **no** es una cadena que se lee igual al revés, sino una cadena igual a su propio **reverso complementario** (A↔T, C↔G). Por ejemplo, `GAATTC` es palíndromo genético, aunque como texto no lo sea. Estas regiones son relevantes porque son propensas a mutaciones.

Por esta definición:
- Solo pueden existir palíndromos de **longitud par**: la base central de uno impar tendría que ser su propio complemento, lo cual es imposible.
- Se adaptó Manacher: se insertan separadores `$` entre las bases, solo se evalúan centros en `$` y la comparación `texto[i - gap] == texto[i + gap]` se reemplaza por una comparación de complementos. La complejidad se mantiene en O(n).

## Archivos de entrada
Se leen automáticamente de la carpeta `data/`:
- `SARS-COV-2-MN908947.3.txt`
- `SARS-COV-2-MT106054.1.txt`
- `gen-M.txt`, `gen-ORF1AB.txt`, `gen-S.txt`
- `seq-proteins.txt`

## Requisitos
- Python 3.x (no requiere librerías externas)

## Uso
Desde la carpeta del proyecto:
```bash
python3 act-integradora1.py
```

## Salidas
En consola:
- Índices y primeros 12 caracteres de cada gen.
- Longitud, posición y secuencia del palíndromo genético más largo por gen.
- Proteínas, índices, primeros 4 aminoácidos y codones.
- Diferencias entre genomas, codones y aminoácidos afectados.
- Cadena más larga sin diferencias.

Archivos que genera el programa en la misma carpeta:
- `palindromos.txt`: palíndromo genético más largo de cada gen (longitud, inicio dentro del gen y secuencia).
- `diferencias.txt`: diferencias entre Wuhan y Texas (índice, base Wuhan, base Texas).

### Resultados de palíndromos
| Gen | Longitud | Inicio dentro del gen | Palíndromo |
|-----|----------|-----------------------|------------|
| M | 10 | 310 | `CGCGTACGCG` |
| S | 16 | 1338 | `GGTAATTATAATTACC` |
| ORF1AB | 20 | 5479 | `ACACTGGTAATTACCAGTGT` |

### Nota sobre los índices
Todos los índices que imprime el programa son **base 0** (Python). NCBI usa base 1, por lo que para comparar contra NCBI hay que sumar 1 (por ejemplo, el gen M aparece en el índice 26522 del programa y en 26523 en NCBI).

## Estructura del repositorio
```
├── act-integradora1.py   # Código fuente
├── data/                 # Secuencias genómicas y proteínas
├── palindromos.txt       # Generado al ejecutar
├── diferencias.txt       # Generado al ejecutar
└── README.md
```

## Entregables
- Programa funcional.
- Video explicativo (~8 min).
- Reporte técnico en equipo.

## Subcompetencias
- SICT0101, SICT0401, STC0101, STC0102.

## Autores
- Andrés Humberto Treviño Garza
- Juan Antonio Rodríguez Reyna
- Juan Esteban Jaramillo Lucero
- Mario Giovanni González López