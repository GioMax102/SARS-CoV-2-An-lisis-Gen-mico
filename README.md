# SARS-CoV-2 Genome Analysis

Análisis computacional del genoma del SARS-CoV-2 mediante algoritmos de procesamiento de strings.

## Descripción
Programa que analiza las secuencias genómicas de Wuhan (2019) y Texas (2020), localiza genes, detecta palíndromos, mapea proteínas y compara diferencias entre ambas versiones del virus.

## Tareas principales
- Localizar índices de los genes **M**, **S** y **ORF1ab** en el genoma.
- Encontrar el palíndromo más largo en cada gen.
- Identificar secciones del genoma que producen cada proteína y sus codones.
- Comparar genomas de Wuhan 2019 vs. Texas 2020: diferencias, codones y aminoácidos afectados.

## Algoritmos
| Punto | Algoritmo | Complejidad |
|-------|-----------|-------------|
| 1. Búsqueda de genes | KMP | O(n + m) |
| 2. Palíndromo más largo | Manacher | O(n) |
| 3. Identificación de proteínas | Traducción en 3 marcos de lectura + KMP | O(n + m) por marco |
| 4. Comparación de genomas | Comparación posición por posición | O(n) |

Todos los algoritmos están implementados desde cero, sin librerías que los resuelvan.

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
- Longitud y secuencia del palíndromo más largo por gen.
- Proteínas, índices, primeros 4 aminoácidos y codones.
- Diferencias entre genomas, codones y aminoácidos afectados.
- Cadena más larga sin diferencias.

Archivos que genera el programa en la misma carpeta:
- `palindromos.txt`: palíndromo más largo de cada gen.
- `diferencias.txt`: diferencias entre Wuhan y Texas (índice, base Wuhan, base Texas).

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