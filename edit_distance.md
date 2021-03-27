# Course 3 in the Bioinformatics Specialization: Comparing Genes, Proteins and Genomes

**Week 2 - Edit Distance**

This function is based on the Smith-Waterman Algorithm and returns the number of edits that has to be made in ```string 1``` for it to become ```string 2```.

```python
def edit_distance(string1, string2):    
    d = {}
    a = len(string1) + 1
    b = len(string2) + 1
    
    for i in range(a):
        d[(i,0)] = i
        
    for j in range(b):
        d[(0,j)] = j
        
    for j in range(1,b):
        for i in range(1,a):
            if string1[i-1] == string2[j-1]:       # if the letters are a match, the score doesn't change
                d[(i,j)] = d[(i-1,j-1)]
            else:                                  # if they are not a match, the score will be increased in 1
                d[(i,j)] = min(d[i-1, j], d[i, j-1], d[i-1, j-1]) + 1
                
    return d[a-1, b-1]                             # returns the score from the last node 
```
    
**Example:**

*Input:*\
string1 = PRTEIN\
string2 = AKIN

*Output:*\
4
