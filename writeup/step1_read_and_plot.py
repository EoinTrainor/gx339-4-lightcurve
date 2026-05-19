"""
Step 1: Read data from a CSV file and plot it as a scatter graph.

Experiment: Photoelectric effect on sodium.
We measured stopping voltage (V_s) for different light frequencies (f).
"""

# --- 1. Import libraries ---
# pandas: reads tables of data (like a spreadsheet)
# matplotlib.pyplot: draws graphs
import pandas as pd
import matplotlib.pyplot as plt

# --- 2. Read the CSV file ---
# pd.read_csv() loads the file into a "DataFrame" — a table with named columns.
# comment='#' tells pandas to skip lines starting with # (our header comments).
data = pd.read_csv("photoelectric_data.csv", comment="#")

# Print the table so we can see what was loaded
print("Loaded data:")
print(data)
print()

# --- 3. Pull out the x and y columns ---
# We access a column by its name in square brackets.
# The column names come from the first non-comment row of the CSV.
x = data["frequency_Hz"]          # x-axis: light frequency in Hz
y = data["stopping_voltage_V"]    # y-axis: stopping voltage in Volts

# --- 4. Create the figure and axes ---
# fig is the whole window; ax is the plot area inside it.
fig, ax = plt.subplots()

# --- 5. Plot the data as a scatter plot ---
# marker='o' draws circles; color and zorder are just styling.
ax.scatter(x, y, marker="o", color="steelblue", zorder=5, label="Measured data")

# --- 6. Label the axes and give the plot a title ---
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Stopping Voltage (V)")
ax.set_title("Photoelectric Effect — Raw Data")

# --- 7. Add a legend and grid ---
ax.legend()
ax.grid(True, linestyle="--", alpha=0.5)

# --- 8. Show the plot ---
plt.tight_layout()   # prevents labels being clipped
plt.show()
