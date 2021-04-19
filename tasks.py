from datetime import datetime
from pathlib import Path
import platform
from invoke import task
from slugify import slugify

SOLUTIONS_DIR = Path("./solutions")
DEFAULT_SOLUTION = "original"
BIN_DIR = "bin"
SRC_DIR = "src"
LOG_DIR = "log"
MACHINE_NAME = slugify(platform.node())
DEFAULT_CC = "gcc"
DEFAULT_CCS = ["gcc-5", "gcc-10", "icc", "clang-6.0", "clang-12", "nvc"]
DEFAULT_CFLAGS = "-O0"
DEFAULT_LDFLAGS = "-lm"
DEFAULT_PARAMFLAGS = "-DN256"
DEFAULT_WFLAGS = "-std=gnu11 -Wall -Wextra -g"
CSV_SEPARATOR = ","
DEFAULT_MEASURE_OPTS = "-e 'cs,migrations,faults' -M 'GFLOPs' "
DEFAULT_RETRIES = 3


def slugify_iter(*args):
    return slugify(" ".join(args))


def get_solution_dir(solution):
    sol_dir = SOLUTIONS_DIR / solution
    if not sol_dir.exists():
        print("Creating dir %s" % sol_dir)
        sol_dir.mkdir(parents=True, exist_ok=True)
    return sol_dir


def get_src_dir(solution):
    src_dir = SOLUTIONS_DIR / solution / SRC_DIR
    if not src_dir.exists():
        print("Creating dir %s" % src_dir)
        src_dir.mkdir(parents=True, exist_ok=True)
    return src_dir


def get_exec_dir(solution):
    exec_dir = SOLUTIONS_DIR / solution / BIN_DIR / MACHINE_NAME
    if not exec_dir.exists():
        print("Creating dir %s" % exec_dir)
        exec_dir.mkdir(parents=True, exist_ok=True)
    return exec_dir


def get_log_dir(solution):
    log_dir = SOLUTIONS_DIR / solution / LOG_DIR / MACHINE_NAME
    if not log_dir.exists():
        print("Creating dir %s" % log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)
    return log_dir


def get_exec_name(cc, cflags, ldflags, paramflags, name="tiny_md"):
    return slugify_iter(name, cc, cflags, ldflags, paramflags)


def get_log_prefix(exec_name, start_time=datetime.now()):
    return "%s-%s" % (start_time.isoformat(timespec="minutes"), exec_name)


def get_flags_str(cc, cflags, ldflags, paramflags, wflags):
    return 'CC="%s" CFLAGS="%s" LDFLAGS="%s" PARAMFLAGS="%s" WFLAGS="%s"' % (
        cc,
        cflags,
        ldflags,
        paramflags,
        wflags,
    )


@task()
def clean_src(c, solution=DEFAULT_SOLUTION):
    """Clean solution source"""
    src_dir = get_src_dir(solution)
    with c.cd(src_dir):
        print(c.cwd)
        print("Cleaning!")
        c.run("make clean")
        print("Cleaned!")


@task()
def clean_exec(c, solution=DEFAULT_SOLUTION):
    """Clean solution executables"""
    exec_dir = get_exec_dir(solution)
    c.run(f"rm -R {exec_dir}")


@task()
def clean_log(c, solution=DEFAULT_SOLUTION):
    """Clean solution exec logs"""
    log_dir = get_log_dir(solution)
    c.run(f"rm -R {log_dir}")


@task()
def clean_solution(c, solution=DEFAULT_SOLUTION):
    """Clean solution"""
    clean_src(c, solution)
    clean_exec(c, solution)
    clean_log(c, solution)


@task(help={"must_clean": "make clean is called prior to build it by default"})
def build(
    c,
    must_clean=True,
    cc=DEFAULT_CC,
    cflags=DEFAULT_CFLAGS,
    wflags=DEFAULT_WFLAGS,
    ldflags=DEFAULT_LDFLAGS,
    paramflags=DEFAULT_PARAMFLAGS,
    solution=DEFAULT_SOLUTION,
):
    """Build project with default options"""
    if must_clean:
        clean_src(c, solution)
    src_dir = get_src_dir(solution)
    with c.cd(src_dir):
        print(c.cwd)
        print("Building!")
        c.run(
            "make tiny_md %s" % get_flags_str(cc, cflags, ldflags, paramflags, wflags)
        )
        print("Built!")
    exec_dir = get_exec_dir(solution)
    exec_path = exec_dir / get_exec_name(cc, cflags, ldflags, paramflags)
    c.run("mv %s/tiny_md %s" % (src_dir, exec_path))
    print("Binary moved to %s" % exec_dir)


@task()
def run(
    c,
    must_clean=True,
    cc=DEFAULT_CC,
    cflags=DEFAULT_CFLAGS,
    wflags=DEFAULT_WFLAGS,
    ldflags=DEFAULT_LDFLAGS,
    paramflags=DEFAULT_PARAMFLAGS,
    solution=DEFAULT_SOLUTION,
    measure_options=DEFAULT_MEASURE_OPTS,
    retries=DEFAULT_RETRIES,
):
    """
    Run a *solution* given some parameters recording a *metric* and
    build it if exec doesn't exist
    """
    exec_dir = get_exec_dir(solution)
    exec_name = get_exec_name(cc, cflags, ldflags, paramflags)
    exec_path = exec_dir / exec_name
    if not exec_path.exists():
        build(c, must_clean, cc, cflags, wflags, ldflags, paramflags, solution)
    log_dir = get_log_dir(solution)
    log_path = log_dir / f"{get_log_prefix(exec_name)}.log"
    stat_csv_path = log_dir / f"{get_log_prefix(exec_name)}-{slugify(measure_options)}.csv"
    command = (
        f"perf stat -x '{CSV_SEPARATOR}' "
        f"{measure_options} "
        f"-r {retries} "
        f"-o {stat_csv_path} "
        f"{exec_path}"
    )
    print("Running $ %s" % command)
    print("Look at log doing $ tail -f %s" % log_path)
    with open(log_path, "a") as log_stream:
        c.run(command, env={"LC_NUMERIC": "en_US"}, out_stream=log_stream)
    c.run(
        f"mv trajectory.xyz {log_dir / get_log_prefix(exec_name)}-trajectory.xyz",
        warn=True,
    )
    c.run(f"mv thermo.log {log_dir / get_log_prefix(exec_name)}-thermo.log", warn=True)


@task()
def run_all_ccs(
    c,
    must_clean=True,
    cflags=DEFAULT_CFLAGS,
    wflags=DEFAULT_WFLAGS,
    ldflags=DEFAULT_LDFLAGS,
    paramflags=DEFAULT_PARAMFLAGS,
    solution=DEFAULT_SOLUTION,
    measure_options=DEFAULT_MEASURE_OPTS,
    retries=DEFAULT_RETRIES,
):
    for cc in DEFAULT_CCS:
        run(
            c,
            must_clean,
            cc,
            cflags,
            wflags,
            ldflags,
            paramflags,
            solution,
            measure_options,
            retries,
        )


@task()
def create_solution(c, name, template=DEFAULT_SOLUTION):
    "Create a new solution with *name* using a *template* src code dir"
    new_src_dir = get_solution_dir(name)
    old_src_dir = get_src_dir(template)
    c.run(f"cp -R {old_src_dir} {new_src_dir}")
    clean_src(c, name)
