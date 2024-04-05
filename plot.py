import csv
import matplotlib.pyplot as plt


csv_file = "./output/energy.csv"


def draw_curve(csv_file):
    x = []
    y = []

    # Reading data from CSV file
    with open(csv_file, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header if present
        for row in reader:
            x.append(float(row[0]))
            y.append(float(row[1]))

    # Plotting the curve
    plt.plot(x, y)
    plt.title('Curve from CSV Data')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    draw_curve(csv_file)
