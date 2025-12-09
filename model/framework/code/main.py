# imports
import os
import csv
import sys
import subprocess
import tempfile
import shutil

# get path of current python executable
python_executable = sys.executable

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

print(smiles_list)

tmp_folder = tempfile.mkdtemp()
tmp_folder = os.path.abspath("tmp")
os.makedirs(tmp_folder, exist_ok=True)
input_tmp = os.path.join(tmp_folder, "input.smi")
with open(input_tmp, "w") as f:
    for i, smi in enumerate(smiles_list):
        f.write(f"{smi} mol_{i}\n")

output_tmp = os.path.join(tmp_folder, "output.csv")

bash_file = os.path.join(tmp_folder, "run_model.sh")
bash_content = f"""
{python_executable} -c "import site; print(site.getsitepackages()[0])"

SITEPKG=$({python_executable} -c "import site; print(site.getsitepackages()[0])")
cat > "$SITEPKG/sitecustomize.py" <<'PY'
try:
    from openbabel import openbabel as ob
    # Provide legacy symbol expected by old code:
    if not hasattr(ob, "OBForceField_FindType") and hasattr(ob, "OBForceField"):
        if hasattr(ob.OBForceField, "FindType"):
            ob.OBForceField_FindType = ob.OBForceField.FindType
        elif hasattr(ob.OBForceField, "FindForceField"):
            ob.OBForceField_FindType = ob.OBForceField.FindForceField
except Exception:
    # Don't break interpreter startup if openbabel isn't installed in some envs
    pass
PY

{python_executable} {root}/entry-cli/calc_props.py -b {input_tmp} -o {output_tmp}
"""

with open(bash_file, "w") as f:
    f.write(bash_content.strip())

# run bash script and print output
process = subprocess.Popen(f"bash {bash_file}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
stdout, stderr = process.communicate()
print("Subprocess output:")
print(stdout)
if stderr:
    print("Subprocess errors:")
    print(stderr)

rb = []
glob = []
primary_amine = []

with open(output_tmp, "r") as f:
    reader = csv.reader(f)
    next(reader)
    for r in reader:
        rb_ = int(r[3])
        glob_ = float(r[4])
        pa_ = str(r[6])

        if rb_ <= 5:
            rb += [1]
        else:
            rb += [0]

        if glob_ <= 0.25:
            glob += [1]
        else:
            glob += [0]

        if pa_ == "True":
            primary_amine += [1]
        else:
            primary_amine += [0]

outputs = []
for i in range(len(smiles_list)):
    r = [rb[i], glob[i], primary_amine[i]]
    outputs += [r]

#check input and output have the same lenght
input_len = len(smiles_list)
output_len = len(outputs)
assert input_len == output_len

# write output in a .csv file
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["rb", "glob", "primary_amine"])
    for o in outputs:
        writer.writerow(o)

# clean up temporary files
shutil.rmtree(tmp_folder)