# Course 3 in the Bioinformatics Specialization: Comparing Genes, Proteins and Genomes

**Week 1 - Manhattan Tourist Problem**

Given a file "matrix.txt" containing the dimensions of the matrix to be analyzed on the first line and two matrices separated by "-" with 
the weights of every edge, such as the one below, this code returns the longest path possible along a matrix with differently weighted edges.

9 10\
3 0 0 4 0 4 1 3 3 0 1\
1 2 2 0 3 2 2 4 4 4 0\
0 2 1 1 0 0 0 3 0 4 2\
2 3 2 0 1 3 2 2 4 3 4\
3 0 4 4 3 2 0 1 1 2 4\
2 3 4 2 1 2 0 2 2 2 2\
4 3 3 4 4 2 2 1 1 2 3\
1 3 0 4 4 3 0 3 2 4 4\
4 1 3 4 0 1 2 3 0 4 2\
*-*\
1 4 1 1 4 3 4 1 0 1\
0 4 4 2 4 3 3 1 2 3\
0 1 2 2 4 2 3 0 2 0\
4 4 4 4 2 3 0 1 1 0\
2 1 3 4 2 0 3 4 1 4\
3 0 2 3 3 3 3 3 0 1\
4 2 1 1 2 1 4 1 2 4\
1 4 0 0 1 4 3 1 3 3\
1 2 3 4 2 0 2 0 4 3\
4 0 4 3 1 2 1 4 1 2

Since we have four informations on this file, we will treat it accordingly.

```python
file = open("matrix.txt").readlines()

lines = [line.rstrip() for line in file]
[n,m] = [int(e) for e in lines[0].split()]  # get the matrix dimensions
b = lines.index('-')
down, right = lines[1:b], lines[b+1:]
down = [[int(e) for e in E.split()] for E in down]  # get the wights of the edges going down
right = [[int(e) for e in E.split()] for E in right]  # get the weight of the edges going right
```

Now we can find the longest possible path in this matrix.

```python
def manhattan_tourist(n, m, down, right):
    s_dic = {(0,0):0}
   
    for i in range(1,n+1):
        s_dic[(i,0)] = s_dic[(i-1,0)] + down[i-1][0]
            
    for j in range(1,m+1):
        s_dic[(0,j)] = s_dic[(0,j-1)] + right[0][j-1]
    
    
    for i in range(1,n+1):
        for j in range(1,m+1):
            s_dic[(i,j)] = max(s_dic[(i-1,j)]+down[i-1][j], s_dic[(i,j-1)] + right[i][j-1])
    
    return s_dic[(n,m)]
 ```
