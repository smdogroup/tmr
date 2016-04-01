# The following code generates the table required for polygonizing a
# face with isosurface values

import matplotlib.pylab as plt
import numpy as np

# This lists out all the cases and indicates the vertex connections
binary = []
vert = []

# None of the corners are active - do nothing
b = [1, 1, 1, 1]
binary.append(b)
vert.append([])

# The face corners - this will use one triangle each
b = [0, 1, 1, 1]
v = [0, 4, 7]
binary.append(b)
vert.append(v)

b = [1, 0, 1, 1]
v = [1, 5, 4]
binary.append(b)
vert.append(v)

b = [1, 1, 0, 1]
v = [2, 6, 5]
binary.append(b)
vert.append(v)

b = [1, 1, 1, 0]
v = [3, 7, 6]
binary.append(b)
vert.append(v)

# Two nodes active - these will use two triangles
b = [0, 0, 1, 1]
v = [0, 1, 5, 0, 5, 7]
binary.append(b)
vert.append(v)

b = [1, 0, 0, 1]
v = [1, 6, 4, 1, 2, 6]
binary.append(b)
vert.append(v)

b = [1, 1, 0, 0]
v = [2, 3, 7, 2, 7, 5]
binary.append(b)
vert.append(v)

b = [0, 1, 1, 0]
v = [0, 4, 3, 3, 4, 6]
binary.append(b)
vert.append(v)

# Corners opposite one another are active
b = [0, 1, 0, 1]
v = [0, 4, 7, 2, 6, 5]
binary.append(b)
vert.append(v)

b = [1, 0, 1, 0]
v = [1, 5, 4, 3, 7, 6] 
binary.append(b)
vert.append(v)

# Three nodes active - these will use 3 triangles
b = [1, 0, 0, 0]
v = [2, 3, 7, 2, 7, 4, 2, 4, 1]
binary.append(b)
vert.append(v)

b = [0, 1, 0, 0]
v = [3, 0, 4, 3, 4, 5, 3, 5, 2]
binary.append(b)
vert.append(v)

b = [0, 0, 1, 0]
v = [0, 1, 5, 0, 5, 6, 0, 6, 3]
binary.append(b)
vert.append(v)

b = [0, 0, 0, 1]
v = [1, 2, 6, 1, 6, 7, 1, 7, 0]
binary.append(b)
vert.append(v)

# All four nodes active
b = [0, 0, 0, 0]
v = [0, 1, 2, 0, 2, 3]
binary.append(b)
vert.append(v)

xpts = np.array([[0, 0], [1, 0], [1, 1], [0, 1],
                 [0.5, 0], [1, 0.5], [0.5, 1], [0, 0.5]])

x = np.array(xpts)
y = 0.0

for k in xrange(len(binary)):
    b = binary[k]
    v = vert[k]

    if k % 4 == 0:
        x = np.array(xpts)
        x[:, 1] += y
        y += 1.25

    for j in xrange(0, len(v), 3):
        plt.plot([x[v[j], 0], x[v[j+1], 0], x[v[j+2], 0], x[v[j], 0]],
                 [x[v[j], 1], x[v[j+1], 1], x[v[j+2], 1], x[v[j], 1]])

    x[:, 0] += 1.25

# Create the vertex and edge tables
edge_table = np.zeros(16, dtype=np.intc)
vert_table = -np.ones((16, 10), dtype=np.intc)

for k in xrange(len(binary)):
    b = binary[k]
    v = vert[k]
    index = 0
    for i in xrange(4):
        if b[i] == 1:
            index |= (1 << i)

    # Index the code
    for i in xrange(len(v)):
        vert_table[index, i] = v[i]
    
    # Compute the vertex code
    code = 0
    for v in vert[k]:
        code |= (1 << v)    

    # Add the code to the edge table
    edge_table[index] = code

# Print out all the results
s = 'const int faceEdgeTable[16] = {%d'%(edge_table[0])
for k in xrange(1, len(binary)):
    s += ', %d'%(edge_table[k])
s += '};\n'

s += 'const int faceTriTable[16][10] = {'
for i in xrange(16):
    s += '\n  {%d'%(vert_table[i, 0])
    for j in xrange(1, 10):
        s += ', %d'%(vert_table[i, j])
    s += '},'
s += '};'

print s

plt.show()
