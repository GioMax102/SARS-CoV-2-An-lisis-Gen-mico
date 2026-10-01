# SARS-CoV-2 Genome Analysis

Análisis computacional del genoma del SARS-CoV-2 mediante algoritmos de procesamiento de strings.

## Descripción
Programa que analiza las secuencias genómicas de Wuhan (2019) y Texas (2020), localiza genes, detecta palíndromos, mapea proteínas y compara diferencias entre ambas versiones del virus.

## Tareas principales
- Localizar índices de los genes **M**, **S** y **ORF1ab** en el genoma.
- Encontrar el palíndromo más largo en cada gen.
- Identificar secciones del genoma que producen cada proteína y sus codones.
- Comparar genomas de Wuhan 2019 vs. Texas 2020: diferencias, codones y aminoácidos afectados.

## Archivos de entrada
- `SARS-COV-2-MN908947.3.txt`
- `SARS-COV-2-MT106054.1.txt`
- `gen-M.txt`, `gen-ORF1AB.txt`, `gen-S.txt`
- `seq-proteins.txt`

## Requisitos
- Python 3.x

## Uso
```bash
python src/main.py
```

## Salidas
- Índices y primeros caracteres de cada gen.
- Longitud del palíndromo más largo por gen.
- Proteínas, índices, primeros aminoácidos y codones.
- Diferencias entre genomas, codones y aminoácidos.
- Cadena más larga sin diferencias.

## Estructura del repositorio
```
├── src/          # Código fuente
├── data/         # Secuencias genómicas
├── output/       # Resultados generados
├── docs/         # Reporte
└── video/        # Presentación
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