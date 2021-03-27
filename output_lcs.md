# Course 3 in the Bioinformatics Specialization: Comparing Genes, Proteins and Genomes

**Week 1 - Longest Commom Subsequence**

This code is divided in two parts. First, the backtrack function. Given two strings `v` and `w`, this function creates a bactrack matrix to store pointers to solve the Longest Common Subsequence Problem:

```python
def lcs_backtrack(v, w):
    s = {(0,0):0}
    backtrack = {}
    
    for i in range(len(v)+1):
        s[(i,0)] = 0               # adds 0 to all Column 0
    
    for j in range(len(w)+1):
        s[(0,j)] = 0               # adds 0 to all Row 0
   
    for a in range(1,len(v)+1):
        for b in range(1,len(w)+1):
            match = 0
            if v[a-1] == w[b-1]:
                match = 1
                
            s[(a,b)] = max(s[(a-1,b)], s[(a,b-1)], s[(a-1,b-1)] + match)      # fills the s matrix with the max possible values based on wheter the letters analyzed are matches
            
            if s[(a,b)] == s[(a-1,b)]:          # creates the backtrack matrix based on the values on the s matrix
                backtrack[(a,b)] = "↓" 
            elif s[(a,b)] == s[(a,b-1)]:
                backtrack[(a,b)] = "→" 
            elif s[(a,b)] == s[(a-1,b-1)] + match:
                backtrack[(a,b)] = "↘" 
                
    return backtrack
 ```
 
 The second part of this problem uses the matrix returned by the `backtrack` function as an argument. The other arguments are the string `v` passed on the `backtrack` function, the length `i` of the string `v` and the length `j` of the string `w`.
 
 ```python
 def output_lcs(backtrack, v, i, j):
    lcs = []
    
    while i > 0 and j > 0:                # when there is a match between the letters of the two strings, this letter is appended to the lcs
        if backtrack[(i,j)] == "↘":
            lcs.append(v[i-1])
            i -= 1
            j -= 1
        elif backtrack[(i,j)] == "↓":
            i -= 1
        elif backtrack[(i,j)] == "→":
            j -= 1
    
    lcs_return = lcs[::-1]                # when appending, the sequence will be backwards, so it is necessary to invert it
    lcs_print = ""
    for letter in lcs_return:
        lcs_print += letter
        
    return lcs_print
 ```

**Examples**

```backtrack``` function\
*Input:*\
PRTEIN\
AKIN

*Output:*\
{(1, 1): '↓',
 (1, 2): '↓',
 (1, 3): '↓',
 (1, 4): '↓',\
 (2, 1): '↓',
 (2, 2): '↓',
 (2, 3): '↓',
 (2, 4): '↓',\
 (3, 1): '↓',
 (3, 2): '↓',
 (3, 3): '↓',
 (3, 4): '↓',\
 (4, 1): '↓',
 (4, 2): '↓',
 (4, 3): '↓',
 (4, 4): '↓',\
 (5, 1): '↓',
 (5, 2): '↓',
 (5, 3): '↘',
 (5, 4): '→',\
 (6, 1): '↓',
 (6, 2): '↓',
 (6, 3): '↓',
 (6, 4): '↘'}
 
 ```output_lcs``` function\
 *Input:*\
 ```lcs_backtrack("PRTEIN", "AKIN"), "PRTEIN", 6, 4```
 
 *Output:*\
 IN
 


 
 
 
  
