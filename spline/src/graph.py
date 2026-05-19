import csv
import matplotlib.pyplot as plt

x_vals = []
f_vals = []
lagrange_vals = []
spline_vals = []

with open("data.csv", "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        x_vals.append(float(row["x"]))
        f_vals.append(float(row["f(x)"]))
        lagrange_vals.append(float(row["lagrange"]))
        spline_vals.append(float(row["spline"]))

xnodes = []
ynodes = []

with open("nodes.csv", "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        xnodes.append(float(row["x"]))
        ynodes.append(float(row["y"]))

plt.figure(figsize=(12, 7))

plt.plot(x_vals, f_vals, linewidth=2, linestyle='-', color="blue", label="исходная функция f(x)")
plt.plot(x_vals, lagrange_vals, linewidth=2,  linestyle='--', color='red', label="интерполяционный многочлен Лагранжа")
plt.plot(x_vals, spline_vals, linewidth=2, linestyle=':', color='green', label="закреплённый кубический сплайн")

plt.scatter(
    xnodes, ynodes,
    s=40,
    marker="^",
    label="узлы интерполяции",
    zorder=5
)

plt.title("интерполяция функции arctg(1/x) - x*x на отрезке [0.1; 1]")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig("plot.png", dpi=200)
plt.show()
