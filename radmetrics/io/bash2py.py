import subprocess

def run_bash_in_python(file_name: str) -> None:
    subprocess.run(['bash', file_name])