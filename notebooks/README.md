# Instrucciones para correr los notebooks

## Prerequisitos

- tener jupyter instalado

## Instalación

1. Crear entorno

        python3 -m venv .env
        source .env/bin/activate
        pip install -r requirements.txt

2. Registrar entorno en jupyter

        python3 -m ipykernel install --user --name=tiny_md_notebooks

3. Iniciar jupyter y elegir el kernel `tiny_md_notebooks` para correr los notebooks
