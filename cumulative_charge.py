import numpy as np
import matplotlib.pyplot as plt
import pandas as pd  # <-- Added for data export

#-----------------------Global Style (Research Paper Look)----------------------
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'legend.fontsize': 12,
    'xtick.direction': 'in',
    'ytick.direction': 'in'
})

#-------------------------Define amino acid categories--------------------------
positive = {'K', 'R'}
negative = {'D', 'E'}
charged = positive.union(negative)

#--------------------------Get category function--------------------------------
def get_category(aa):
    if aa in positive:
        return 'positive'
    elif aa in negative:
        return 'negative'
    else:
        return 'neutral'

#--------------------------------Input FASTA-----------------------------------------------------------
print("Paste your FASTA sequence (starting with >header is optional). Press Enter twice when done:\n")
fasta_lines = []
while True:
    try:
        line = input()
        if not line.strip():
            break
        fasta_lines.append(line.strip())
    except EOFError:
        break

#-------------------------------Parse FASTA-------------------------------------
if fasta_lines and fasta_lines[0].startswith(">"):
    sequence_name = fasta_lines[0][1:].strip()
    sequence = ''.join(fasta_lines[1:]).replace(" ", "").upper()
else:
    sequence_name = "nup100"
    sequence = ''.join(fasta_lines).replace(" ", "").upper()

#--------Count and accumulate positive, negative, and charged residues----------
positions = []
charged_step_curve = []
positive_step_curve = []
negative_step_curve = []

charged_count = 0
positive_count = 0
negative_count = 0

for i, aa in enumerate(sequence):
    if aa in charged:
        charged_count += 1
    if aa in positive:
        positive_count += 1
    if aa in negative:
        negative_count += 1

    positions.append(i + 1)
    charged_step_curve.append(charged_count)
    positive_step_curve.append(positive_count)
    negative_step_curve.append(negative_count)

#---------------------------------Data export for OriginPro---------------------
data = {
    "Residue_Position": positions,
    "Charged(D,E,K,R)": charged_step_curve,
    "Positive(K,R)": positive_step_curve,
    "Negative(D,E)": negative_step_curve
}
df = pd.DataFrame(data)

data_file = f"{sequence_name}_charge__data.csv"
df.to_csv(data_file, index=False)
print(f"Data file for OriginPro saved as: {data_file}")

#-----------------------------------------Plotting------------------------------
fig, ax = plt.subplots(figsize=(10, 6))

#----------------------------------Charged (black, thick line)------------------
ax.plot(positions, charged_step_curve, label='Charged (D, E, K, R)',
        color='black', linewidth=3)

#-------------------------------------Positive (blue)---------------------------
ax.plot(positions, positive_step_curve, label='Positive (K, R)',
        color='blue', linewidth=2)

#---------------------------------------Negative (red)--------------------------
ax.plot(positions, negative_step_curve, label='Negative (D, E)',
        color='red', linewidth=2)

#---------------Optional vertical line (rigid vs flexible boundary)-------------
ax.axvline(x=610, color='gray', linestyle='--', linewidth=1.5)

#----------------------Labels and aesthetics------------------------------------
ax.set_xlabel('Residue Position')
ax.set_ylabel('Cumulative Count')
#ax.set_title(f'Cumulative Charge Residues Distribution along sequences:NUP100')
ax.legend(frameon=False, loc='upper left')

ax.grid(True, linestyle='--', linewidth=0.6, alpha=0.7)
ax.set_facecolor('white')

plt.tight_layout()

#-----------------------------------Save High-Quality Figures-------------------
# pdf_file = f"{sequence_name}_charge_distribution.pdf"
png_file = f"{sequence_name}_curve1.png"

plt.savefig(png_file, dpi=600, bbox_inches='tight')  # High-res raster
# plt.savefig(pdf_file, dpi=600, bbox_inches='tight')   # Uncomment if needed

#----------------------------Show Plot------------------------------------------
plt.show()

print(f"Saved high-quality figure as:\n- {png_file}")

