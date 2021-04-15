# tiny_md

> _tiny molecular dynamics_

Proyecto inicial para realizar speed-up de un código de Dinámica Molecular con un potencial interatómico de Lennard-Jones. Inicialmente se dan las posiciones en un cristal FCC y velocidades aleatorias distribuidas uniformemente según la temperatura inicial. Las fuerzas se obtienen de un potencial de LJ (12-6) y la evolución temporal viene dada por el algoritmo Velocity Verlet. Se consideran condiciones periódicas de contorno para reproducir un sistema infinito. La temperatura es controlada a través de un reescaleo en las velocidades. Cada cierta cantidad de pasos de dinámica molecular se cambia la densidad del sistema y se reescalean las posiciones para obtener la ecuación de estado. Para más información se puede ver `TinyMD.pdf`.


### Requisitos

Para compilar es necesario tener instalado `gcc` y `OpenGL`.


### Compilación

Para compilar se utiliza `Makefile`:
```bash
make clean
make
```
donde `make clean` elimina los objetos compilados anteriormente y `make` compila dos ejecutables: `tiny_md` y `viz`, ambos realizan la misma simulación pero el segundo posee una visualización en tiempo real.

> Nota:
>
> _Si se desean cambiar parámetros de entrada de la simulación, puede modificarse el archivo _`parameters.h`_ o pasar los valores deseados como parámetros al preprocesador C; por ejemplo, _`make CPPFLAGS="-DN=1372"`_ cambia la cantidad de partículas que se simulan._

### Archivos
*Quizás darle otra ubicación a esta lista*
* core.c : contiene todas las funciones: calcula posic. y veloc. inic., imagen mínima, fuerzas, PBC, y pasos con Velocity Verlet.
* core.h : declara funciones usadas por main, core.c, etc.
* Makefile : para compilar.
* meson.build : bien gracias.
* parameters.h : contiene parámetros: nro. partículas, T, densidad, radio de corte, etc.
* README.md : este archivo.
* tiny\_md.c : código princ. de la DM.
* TinyMD.pdf : sobre DM y algoritmos usados en el programa.
* viz.c : como `tiny\_md.c` con visualización simultánea.
* wtime.c : revisar: sobre calcular t sim.
* wtime.h : `#pragma once`
### Contacto

Por errores, preguntas o sugerencias contactarse con:

+ Francisco Fernandez (<fernandezfrancisco2195@gmail.com>)
