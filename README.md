# tiny_md

> _tiny molecular dynamics_

Proyecto inicial para realizar speed-up de un código de Dinámica Molecular con un potencial interatómico de Lennard-Jones. Inicialmente se dan las posiciones en un cristal FCC y velocidades aleatorias distribuidas uniformemente según la temperatura inicial. Las fuerzas se obtienen de un potencial de LJ (12-6) y la evolución temporal viene dada por el algoritmo Velocity Verlet. Se consideran condiciones periódicas de contorno para reproducir un sistema infinito. La temperatura es controlada a través de un reescaleo en las velocidades. Cada cierta cantidad de pasos de dinámica molecular se cambia la densidad del sistema y se reescalean las posiciones para obtener la ecuación de estado. Para más información se puede ver `TinyMD.pdf`.

## Prerequisitos

* python>=3.6

    `sudo apt install python3 python3-pip` en Ubuntu 18.04

* invoke==1.5.0
  * recomendación instalar con [pipx](https://pipxproject.github.io/pipx/installation/) como user (sin `sudo`), el paquete `python3-invoke` suele estar desactualizado

    ```bash
    python3 -m pip install --user pipx
    python3 -m pipx ensurepath
    pipx install invoke==1.5.0
    pipx inject invoke python-slugify
    ```
  
* perf
  * `sudo apt install --install-recommends linux-generic-hwe-18.04 xserver-xorg-hwe-18.04 linux-tools-generic-hwe-18.04` en Ubuntu 18.04 para usar el último kernel.
  
    Sino `sudo apt install linux-tools-generic`
  * También es necesario que `cat /proc/sys/kernel/perf_event_paranoid` sea -1.

    Si no es -1, hacer `echo -1 | sudo tee /proc/sys/kernel/perf_event_paranoid`

### Hyperthreading

Verificar si Hyperthreading funciona

    cat /sys/devices/system/cpu/smt/active

Apagar hyperthreading

    echo off | sudo tee /sys/devices/system/cpu/smt/control

Volver a prender

    echo on | sudo tee /sys/devices/system/cpu/smt/control


## Como correr soluciones

El proyecto esta estructurado en soluciones (`solutions` dir). Cada solucion tiene su carpeta y dentro de esta su fuente (`src`), sus binarios executables (`bin`) y sus registros (`log`), estos 2 últimos son guardados a su vez en una carpeta por máquina

    solutions/$name/src
    solutions/$name/bin/$machine/
    solutions/$name/log/$machine/

El código original está en la solución `original`

Para crear una nueva solución basada en la original

    invoke create-solution $name

Nota: se puede agregar `--template=$existint_solution_name` para copiar de otra solucion que no sea la original

Para ejecutar una solucion

    invoke run --solution $name

Esto compila creando un nuevo ejecutable y lo ejecuta.

Para ver todas las opciones de ejecucion mirar:

    invoke run --help

Para ver todos los comandos hacer

    invoke --list

## Requisitos

Para compilar es necesario tener instalado:

* `gcc`
* `OpenGL`
  * Instalar con `sudo apt install freeglut3-dev` en Ubuntu 18.04

## Compilación

Para compilar se utiliza `Makefile`:

```bash
make clean
make
```

donde `make clean` elimina los objetos compilados anteriormente y `make` compila dos ejecutables: `tiny_md` y `viz`, ambos realizan la misma simulación pero el segundo posee una visualización en tiempo real.

> Nota:
>
> _Si se desean cambiar parámetros de entrada de la simulación, puede modificarse el archivo _`parameters.h`_ o pasar los valores deseados como parámetros al preprocesador C; por ejemplo, _`make PARAMFLAGS="-DN=1372"`_ cambia la cantidad de partículas que se simulan._

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

Por errores, preguntas o sugerencias contactarse con autor original :D

+ Francisco Fernandez (<fernandezfrancisco2195@gmail.com>)
