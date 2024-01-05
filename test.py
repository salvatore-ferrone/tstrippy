from pathlib import Path
this_directory = Path(__file__).parent
print(this_directory)
long_description = (this_directory / "README.md").read_text()
print(long_description)