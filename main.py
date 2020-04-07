import random, math, numpy, copy
from scipy.stats import t, f

N = 8
m = 3
p = 0.85
q = 1-p
d = 4
f1 = m-1
f2 = N
f3 = f1*f2
f4 = N-d


x1min = 10
x1max = 60
x2min = -70
x2max = -10
x3min = 60
x3max = 70
xAvmin = (x1min + x2min + x3min)/3
xAvmax = (x1max + x2max + x3max)/3
ymin = 200 + xAvmin
ymax = 200 + xAvmax
random.seed()

table_NormExperiment = [["N", "x0", "x1", "x2", "x3", "x1x2", "x1x3", "x2x3", "x1x2x3"],
                        [1,     1,   -1,   -1,   -1,     1,      1,      1,      -1],
                        [2,     1,   -1,   -1,    1,     1,     -1,     -1,       1],
                        [3,     1,   -1,    1,   -1,    -1,      1,     -1,       1],
                        [4,     1,   -1,    1,    1,    -1,     -1,      1,      -1],
                        [5,     1,    1,   -1,   -1,    -1,     -1,      1,       1],
                        [6,     1,    1,   -1,    1,    -1,      1,     -1,      -1],
                        [7,     1,    1,    1,   -1,     1,     -1,     -1,      -1],
                        [8,     1,    1,    1,    1,     1,      1,      1,       1]]

table_NaturExperiment = [["N", "x0", "x1", "x2", "x3", "x1x2",       "x1x3",       "x2x3",        "x1x2x3"],
                         [1,     1, x1min, x2min, x3min, x1min*x2min,  x1min*x3min,  x2min*x3min,  x1min*x2min*x3min],
                         [2,     1, x1min, x2min, x3max, x1min*x2min,  x1min*x3max,  x2min*x3max,  x1min*x2min*x3max],
                         [3,     1, x1min, x2max, x3min, x1min*x2max,  x1min*x3min,  x2max*x3min,  x1min*x2max*x3min],
                         [4,     1, x1min, x2max, x3max, x1min*x2max,  x1min*x3max,  x2max*x3max,  x1min*x2max*x3max],
                         [5,     1, x1max, x2min, x3min, x1max*x2min,  x1max*x3min,  x2min*x3min,  x1max*x2min*x3min],
                         [6,     1, x1max, x2min, x3max, x1max*x2min,  x1max*x3max,  x2min*x3max,  x1max*x2min*x3max],
                         [7,     1, x1max, x2max, x3min, x1max*x2max,  x1max*x3min,  x2max*x3min,  x1max*x2max*x3min],
                         [8,     1, x1max, x2max, x3max, x1max*x2max,  x1max*x3max,  x2max*x3max,  x1max*x2max*x3max]]

coef_eqFull = [[1, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0]]

coef_eqLin = [[1, 0, 0, 0],
              [0, 0, 0, 0],
              [0, 0, 0, 0],
              [0, 0, 0, 0]]


free_el = [0, 0, 0, 0, 0, 0, 0, 0]

for i in range(8):
    for j in range(max(i, 1), 8):
        for n in range(1, 9):
            coef_eqFull[i][j] += table_NaturExperiment[n][i+1]*table_NaturExperiment[n][j+1]
        coef_eqFull[i][j] /= 8
        coef_eqFull[j][i] = coef_eqFull[i][j]
        if i < 4 and j < 4:
            coef_eqLin[i][j] = coef_eqFull[i][j]
            coef_eqLin[j][i] = coef_eqFull[j][i]


b = [0, 0, 0, 0, 0, 0, 0, 0]
a = [0, 0, 0, 0, 0, 0, 0, 0]

AvYs = []

DisYs = []

t_val = []
t_cr = 0.0
insign = []

F_val = 0.0
F_cr = 0.0

def print_table(table):
    if not "Yi_av" in table_NormExperiment[0]:
        table_NormExperiment[0].append("Yi_av")
        table_NormExperiment[0].append("S_Yi")
    print("\n", "-" * len(table_NormExperiment[0]) * 11, "\n")
    for i in range(9):
        print("|", end="")
        for j in range(len(table_NormExperiment[i])):
            if i > 0 and j > 8:
                print("{:.2f}".format(float(table_NormExperiment[i][j])), end="    |")
            elif i == 0 and j > 8:
                print(table_NormExperiment[i][j], " "*(9-len(str(table_NormExperiment[i][j]))), end="|")
            else:
                print(table[i][j], " " * (9 - len(str(table[i][j]))), end="|")
        if i > 0:
            print("{:.2f}".format(float(AvYs[i-1])), end="    |")
            print("{:.2f}".format(float(DisYs[i-1])), end="    |")
        print("\n", "-" * len(table_NormExperiment[0])*11, "\n")


def randomize(s, e):
    global AvYs
    AvYs = []
    global DisYs
    DisYs = []
    for i in range(9):
        sum_y = 0
        for j in range(s, e):
            if i == 0:
                table_NormExperiment[i].append("Yi{}".format(j - 2))
            else:
                y = random.uniform(ymin, ymax)
                table_NormExperiment[i].append(y)
                sum_y += y
        if not i == 0:
            AvYs.append(sum_y/m)
    for i in range(1, 9):
        sum_y = 0
        for j in range(3, m+3):
            sum_y += pow(table_NormExperiment[i][j]- AvYs[i-1], 2)
        DisYs.append(sum_y/(m-1))


def cochran():
    global DisYs
    max_dispersion = max(DisYs)
    Gp = max_dispersion/sum(DisYs)
    fisher = table_fisher(p, 1, f3)
    Gt = fisher/(fisher+f2-1)
    return Gp < Gt


def table_fisher(prob, d, f3):
    x_vec = [i*0.001 for i in range(int(10/0.001))]
    for i in x_vec:
        if abs(f.cdf(i, N-d, f3)-prob) < 0.0001:
            return i


def coef(lin):
    global a
    global b
    b = [0, 0, 0, 0, 0, 0, 0, 0]
    a = [0, 0, 0, 0, 0, 0, 0, 0]
    free_el[0] = sum(AvYs)/len(AvYs)
    b[0] = free_el[0]
    for i in range(1, 8):
        for n in range(1, 9):
            free_el[i] += table_NaturExperiment[n][i+1]*AvYs[n-1]
            b[i] += table_NormExperiment[n][i+1]*AvYs[n-1]
        free_el[i] /= 8
        b[i] /= 8
    if lin:
        denominator = numpy.linalg.det(numpy.array(coef_eqLin))
        for i in range(4):
            numerM = copy.deepcopy(coef_eqLin)
            for j in range(4):
                numerM[j][i] = free_el[j]
            numerator = numpy.linalg.det(numpy.array(numerM))
            a[i] = numerator/denominator

    else:
        denominator = numpy.linalg.det(numpy.array(coef_eqFull))
        for i in range(8):
            numerM = copy.deepcopy(coef_eqFull)
            for j in range(8):
                numerM[i][j] = free_el[j]
            numerator = numpy.linalg.det(numpy.array(numerM))
            a[i] = numerator / denominator


def student(lin):
    global DisYs
    global AvYs
    global d
    global t_val
    global t_cr
    global insign
    global b
    AvDisYs = sum(DisYs)/len(DisYs)
    Sb = math.sqrt(AvDisYs/(N*m))
    t_val = []
    if lin:
        r = 4
    else:
        r = 8
    for x in range(r):
        t_val.append(math.fabs(b[x])/Sb)
    t_cr = 0
    x_vec = [i * 0.0001 for i in range(int(5 / 0.0001))]
    par = 0.5 + p / 0.1 * 0.05
    for i in x_vec:
        if abs(t.cdf(i, f3) - par) < 0.000005:
            t_cr = i
            break
    insign = []
    for i in range(len(t_val)):
        if t_val[i] <= t_cr:
            insign.append(i)
            d -= 1


def fisher():
    global DisYs
    global F_val
    global F_cr
    AvDisYs = sum(DisYs) / len(DisYs)
    Sad = 0
    for dis in DisYs:
        Sad += dis*(m-1)
    Sad = Sad*m/(N-d)
    F_val = Sad/AvDisYs
    x_vec = [i * 0.001 for i in range(int(10 / 0.001))]
    F_cr = None
    for i in x_vec:
        if abs(f.cdf(i, N - d, f3) - p) < 0.0001:
            F_cr = i
            break
    if not F_cr:
        print("\nSomething went wrong.\nUnable to calculate critical value for Fisher's test")
    elif F_cr >= F_val:
        print("\nF = {}\t\t\tF_cr = {}\t\t\tF =< F_cr\nAccording to Fisher's F-test model is adequate to the original.".format(F_val, F_cr))
        return True
    else:
        print("\nF = {}\t\t\tF_cr = {}\t\t\tF > F_cr\nAccording to Fisher's F-test model is not adequate to the original.".format(F_val, F_cr))
        return False


def printRes(lin):
    print("Normalized Experiment:")
    print_table(table_NormExperiment)
    print("\nNaturalized Experiment:")
    print_table(table_NaturExperiment)
    print("\n\nNormalized equation:\ny = ", end="")
    for i in range(4):
        if i == 0:
            print("{:.2f} ".format(b[i]), end="")
        else:
            print("{:+.2f}*x{}".format(b[i], i), end="")
    if not lin:
        print("{:+.2f}*x1x2 {:+.2f}*x1x3 {:+.2f}*x2x3 {:+.2f}*x1x2x3".format(b[4], b[5], b[6], b[7]))
        print("\nCheck:")
        for i in range(1, 8):
            print("{:.2f} {:+.2f} {:+.2f} {:+.2f} {:+.2f} {:+.2f} {:+.2f} {:+.2f} = {:.2f}\ny = {:.2f}\n".format(b[0],
                                                                                 b[1] * table_NormExperiment[i][2],
                                                                                 b[2] * table_NormExperiment[i][3],
                                                                                 b[3] * table_NormExperiment[i][4],
                                                                                 b[4] * table_NormExperiment[i][5],
                                                                                 b[5] * table_NormExperiment[i][6],
                                                                                 b[6] * table_NormExperiment[i][7],
                                                                                 b[7] * table_NormExperiment[i][8],
                                                                                 b[0] +
                                                                                 b[1] * table_NormExperiment[i][2] +
                                                                                 b[2] * table_NormExperiment[i][3] +
                                                                                 b[3] * table_NormExperiment[i][4] +
                                                                                 b[4] * table_NormExperiment[i][5] +
                                                                                 b[5] * table_NormExperiment[i][6] +
                                                                                 b[6] * table_NormExperiment[i][7] +
                                                                                 b[7] * table_NormExperiment[i][8],
                                                                                 AvYs[i - 1]))
    else:
        print("\nCheck:")
        for i in range(1, 8):
            print("{:.2f} {:+.2f} {:+.2f} {:+.2f} = {:.2f}\ny = {:.2f}\n".format(b[0],
                                                                                 b[1] * table_NormExperiment[i][2],
                                                                                 b[2] * table_NormExperiment[i][3],
                                                                                 b[3] * table_NormExperiment[i][4],
                                                                                 b[0] +
                                                                                 b[1] * table_NormExperiment[i][2] +
                                                                                 b[2] * table_NormExperiment[i][3] +
                                                                                 b[3] * table_NormExperiment[i][4],
                                                                                 AvYs[i - 1]))
    print("\n\nNaturalized equasion:\ny = ", end="")
    for i in range(4):
        if i == 0:
            print("{:.2f} ".format(a[i]), end="")
        else:
            print("{:+.2f}*x{}".format(a[i], i), end="")
    if not lin:
        print("{:+.2f}*x1x2 {:+.2f}*x1x3 {:+.2f}*x2x3 {:+.2f}*x1x2x3".format(a[4], a[5], a[6], a[7]))
        print("\nCheck:")
        for i in range(1, 8):
            print("{:.2f} {:+.2f} {:+.2f} {:+.2f} {:+.2f} {:+.2f} {:+.2f} {:+.2f} = {:.2f}\ny = {:.2f}\n".format(b[0],
                                                                                 a[1] * table_NaturExperiment[i][2],
                                                                                 a[2] * table_NaturExperiment[i][3],
                                                                                 a[3] * table_NaturExperiment[i][4],
                                                                                 a[4] * table_NaturExperiment[i][5],
                                                                                 a[5] * table_NaturExperiment[i][6],
                                                                                 a[6] * table_NaturExperiment[i][7],
                                                                                 a[7] * table_NaturExperiment[i][8],
                                                                                 a[0] +
                                                                                 a[1] * table_NaturExperiment[i][2] +
                                                                                 a[2] * table_NaturExperiment[i][3] +
                                                                                 a[3] * table_NaturExperiment[i][4] +
                                                                                 a[4] * table_NaturExperiment[i][5] +
                                                                                 a[5] * table_NaturExperiment[i][6] +
                                                                                 a[6] * table_NaturExperiment[i][7] +
                                                                                 a[7] * table_NaturExperiment[i][8],
                                                                                 AvYs[i - 1]))
    else:
        print("\nCheck:")
        for i in range(1, 8):
            print("{:.2f} {:+.2f} {:+.2f} {:+.2f} = {:.2f}\ny = {:.2f}\n".format(b[0],
                                                                                 a[1] * table_NaturExperiment[i][2],
                                                                                 a[2] * table_NaturExperiment[i][3],
                                                                                 a[3] * table_NaturExperiment[i][4],
                                                                                 a[0] +
                                                                                 a[1] * table_NaturExperiment[i][2] +
                                                                                 a[2] * table_NaturExperiment[i][3] +
                                                                                 a[3] * table_NaturExperiment[i][4],
                                                                                 AvYs[i - 1]))
    print("According to Student's t-test these coefficients are insignificant:")
    for ind in insign:
        print("t = {}\t\t\tt_cr = {}\t\t\tt < t_cr\nb{} = {:.2f} and a{} = {:.2f}".format(t_val[ind], t_cr, ind, b[ind], ind,a[ind]))
    print("\nThen the equations change:\nNormalized:\ny = ", end="")
    for i in range(4):
        if not i in insign:
            if i == 0:
                print("{:.2f} ".format(b[i]), end="")
            else:
                print("{:+.2f}*x{}".format(b[i], i), end="")
    if not lin:
        for i in range(4, 8):
            if not i in insign:
                print("{:+.2f}*{}".format(b[i], table_NormExperiment[0][i+1]), end="")
    print("\n\nNaturalized:\ny = ", end="")
    for i in range(4):
        if not i in insign:
            if i == 0:
                print("{:.2f} ".format(a[i]), end="")
            else:
                print("{:+.2f}*x{}".format(a[i], i), end="")
    if not lin:
        for i in range(4, 8):
            if not i in insign:
                print("{:+.2f}*{}".format(a[i], table_NormExperiment[0][i+1]), end="")
    if not F_cr:
        print("\n\nSomething went wrong.\nUnable to calculate critical value for Fisher's test")
    elif F_cr >= F_val:
        print("\n\nF = {}\t\t\tF_cr = {}\t\t\tF =< F_cr\nAccording to Fisher's F-test model is adequate to the original.".format(F_val, F_cr))
    else:
        print("\n\nF = {}\t\t\tF_cr = {}\t\t\tF > F_cr\nAccording to Fisher's F-test model is not adequate to the original.".format(F_val, F_cr))





exp_counter = 0
fish = False

while exp_counter < 10 and not fish:
    print("\nTest set #{}\t\tLinear equation\nGenerating values...".format(exp_counter+1))
    m = 3
    N = 4
    d = 4
    f1 = m - 1
    f4 = N - d
    startY = 3
    endY = m + 3
    table_NormExperiment = [["N", "x0", "x1", "x2", "x3", "x1x2", "x1x3", "x2x3", "x1x2x3"],
                            [1, 1, -1, -1, -1, 1, 1, 1, -1],
                            [2, 1, -1, -1, 1, 1, -1, -1, 1],
                            [3, 1, -1, 1, -1, -1, 1, -1, 1],
                            [4, 1, -1, 1, 1, -1, -1, 1, -1],
                            [5, 1, 1, -1, -1, -1, -1, 1, 1],
                            [6, 1, 1, -1, 1, -1, 1, -1, -1],
                            [7, 1, 1, 1, -1, 1, -1, -1, -1],
                            [8, 1, 1, 1, 1, 1, 1, 1, 1]]

    table_NaturExperiment = [["N", "x0", "x1", "x2", "x3", "x1x2", "x1x3", "x2x3", "x1x2x3"],
                             [1, 1, x1min, x2min, x3min, x1min * x2min, x1min * x3min, x2min * x3min,
                              x1min * x2min * x3min],
                             [2, 1, x1min, x2min, x3max, x1min * x2min, x1min * x3max, x2min * x3max,
                              x1min * x2min * x3max],
                             [3, 1, x1min, x2max, x3min, x1min * x2max, x1min * x3min, x2max * x3min,
                              x1min * x2max * x3min],
                             [4, 1, x1min, x2max, x3max, x1min * x2max, x1min * x3max, x2max * x3max,
                              x1min * x2max * x3max],
                             [5, 1, x1max, x2min, x3min, x1max * x2min, x1max * x3min, x2min * x3min,
                              x1max * x2min * x3min],
                             [6, 1, x1max, x2min, x3max, x1max * x2min, x1max * x3max, x2min * x3max,
                              x1max * x2min * x3max],
                             [7, 1, x1max, x2max, x3min, x1max * x2max, x1max * x3min, x2max * x3min,
                              x1max * x2max * x3min],
                             [8, 1, x1max, x2max, x3max, x1max * x2max, x1max * x3max, x2max * x3max,
                              x1max * x2max * x3max]]
    randomize(startY, endY)
    eqlin = True
    cochran_cond = cochran()
    while not cochran_cond:
        m += 1
        startY = endY
        endY = m + 3
        f1 = m - 1
        randomize(startY, endY)
        cochran_cond = cochran()
    print("Calculating coefficients...")
    coef(eqlin)
    student(eqlin)
    f4 = N - d
    fish = fisher()
    if not fish:
        print("\t\t\tNon-linear equation\nGenerating values...")
        m = 3
        N = 8
        d = 8
        f1 = m - 1
        f4 = N - d
        startY = 3
        endY = m + 3
        eqlin = False
        table_NormExperiment = [["N", "x0", "x1", "x2", "x3", "x1x2", "x1x3", "x2x3", "x1x2x3"],
                                [1, 1, -1, -1, -1, 1, 1, 1, -1],
                                [2, 1, -1, -1, 1, 1, -1, -1, 1],
                                [3, 1, -1, 1, -1, -1, 1, -1, 1],
                                [4, 1, -1, 1, 1, -1, -1, 1, -1],
                                [5, 1, 1, -1, -1, -1, -1, 1, 1],
                                [6, 1, 1, -1, 1, -1, 1, -1, -1],
                                [7, 1, 1, 1, -1, 1, -1, -1, -1],
                                [8, 1, 1, 1, 1, 1, 1, 1, 1]]

        table_NaturExperiment = [["N", "x0", "x1", "x2", "x3", "x1x2", "x1x3", "x2x3", "x1x2x3"],
                                 [1, 1, x1min, x2min, x3min, x1min * x2min, x1min * x3min, x2min * x3min,
                                  x1min * x2min * x3min],
                                 [2, 1, x1min, x2min, x3max, x1min * x2min, x1min * x3max, x2min * x3max,
                                  x1min * x2min * x3max],
                                 [3, 1, x1min, x2max, x3min, x1min * x2max, x1min * x3min, x2max * x3min,
                                  x1min * x2max * x3min],
                                 [4, 1, x1min, x2max, x3max, x1min * x2max, x1min * x3max, x2max * x3max,
                                  x1min * x2max * x3max],
                                 [5, 1, x1max, x2min, x3min, x1max * x2min, x1max * x3min, x2min * x3min,
                                  x1max * x2min * x3min],
                                 [6, 1, x1max, x2min, x3max, x1max * x2min, x1max * x3max, x2min * x3max,
                                  x1max * x2min * x3max],
                                 [7, 1, x1max, x2max, x3min, x1max * x2max, x1max * x3min, x2max * x3min,
                                  x1max * x2max * x3min],
                                 [8, 1, x1max, x2max, x3max, x1max * x2max, x1max * x3max, x2max * x3max,
                                  x1max * x2max * x3max]]
        randomize(startY, endY)
        cochran_cond = cochran()
        while not cochran_cond:
            m += 1
            startY = endY
            endY = m + 3
            f1 = m - 1
            randomize(startY, endY)
            cochran_cond = cochran()
        print("Calculating coefficients...")
        coef(eqlin)
        student(eqlin)
        f4 = N - d
        fish = fisher()
    exp_counter += 1

if fish:
    printRes(eqlin)
else:
    print("\nNeither of 10 sets passed the Fisher test.\nThe results of the last test set:\n")
    printRes(eqlin)
