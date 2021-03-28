# Course 3 in the Bioinformatics Specialization: Comparing Genes, Proteins and Genomes

**Week 2 - Overlap Alignment**

This function finds a suffix on ``` string 1``` that best overlaps with a prefix on ```string 2```, creating an alignment that overlaps the two strings.

```python
def overlap_alignment(string1, string2, premium, penalty):
    score_matrix = {}
    backtrack = {}
    
    for i in range(len(string1)+1):
        score_matrix[(i,0)] = 0                # adds 0 to Column 0
        
    for j in range(len(string2)+1):
        score_matrix[(0,j)] = j * penalty      # adds j * penalty to Row 0
    
    location = 0
    for index, letter in enumerate(string1):   # finds the firt time the first letter from string 2 appears in string 1, this will be where the alignment will begin
        if letter == string2[0]:
            location = index
            break
    
    for i in range(location+1):                 # add 0 to all rows in the first column until the row where the first letter in string 2 appears for the first time in string 1
        score_matrix[(i,1)] = 0
        
    for j in range(len(string2)+1):             # add 0 to all columns in the rows preceding the row where the first letter in string 2 appear for the first time in string 1
      for i in range(1,location+1):
            score_matrix[(i,j)] = 0
    
    for a in range(location+1, len(string1)+1):         # adds to the score matrix a premium value if match and penalty value if mismatch or indel
        for b in range(1, len(string2)+1):
            match = premium * -10
            mismatch = penalty * 10
            if string1[a-1] == string2[b-1]:
                match = score_matrix[(a-1, b-1)] + premium
            elif string1[a-1] != string2[b-1]:
                mismatch = score_matrix[(a-1, b-1)] + penalty
            delete = score_matrix[(a-1,b)] + penalty
            insert = score_matrix[(a,b-1)] + penalty
            score_matrix[(a,b)] = max(match, mismatch, insert, delete)       # creates backtrack matrix based on highest scores from score_matrix
            if score_matrix[(a,b)] == delete:
                backtrack[(a,b)] = "d"
            elif score_matrix[(a,b)] == insert:
                backtrack[(a,b)] = "r"
            elif score_matrix[(a,b)] == mismatch:
                backtrack[(a,b)] = "m"
            elif score_matrix[(a,b)] == match:
                backtrack[(a,b)] = "c"
                
    max_node = 0
    max_node_score = 0
    
    for j in range(len(string2)+1):             # find first node in last row where value is max to use it as starting point for backtrakcing 
        if score_matrix[(len(string1),j)] > max_node_score:
            max_node = j
            max_node_score = score_matrix[(len(string1),j)] 
    
    i = len(string1)
    j = max_node
    alig_a = ""
    alig_b = ""
    final_score = 0
    
    while i > location and j > 0:               # backtracking and construction of alignments
        
        if backtrack[(i,j)] == "c":
            alig_a = string1[i-1] + alig_a
            alig_b = string2[j-1] + alig_b
            i -= 1
            j -= 1
            final_score += premium
            
        elif backtrack[(i,j)] == "r":
            alig_a = "-" + alig_a
            alig_b = string2[j-1] + alig_b
            j -= 1
            final_score += penalty
            
        elif backtrack[(i,j)] == "d":
            alig_a = string1[i-1] + alig_a
            alig_b = "-" + alig_b
            i -= 1
            final_score += penalty
        else:
            alig_a = string1[i-1] + alig_a
            alig_b = string2[j-1] + alig_b
            i -= 1
            j -= 1
            final_score += penalty
        
    return final_score, alig_a, alig_b
```

**Example**\
*Input:*\
PRTEINS\
FEING
1\
-1

*Output:*\
 1, '-EINS', 'FEIN-'  
    
    
    
    
