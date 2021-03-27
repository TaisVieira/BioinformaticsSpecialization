# Course 3 in the Bioinformatics Specialization: Comparing Genes, Proteins and Genomes

**Week 2 - Global Alignment**

This code is based on the Needleman Wunsch algorithm. Usign a `scoring matrix` for matches/mismatches and a `penalty` variable for insertions or deletions, the function returns the highest possible score for a **global** alignment between the two strings along with what this alignment is.

```python
def global_alignment(string1, string2, penalty, score):
    d = penalty
    score_matrix = {}
    
    for i in range(len(string1)+1):
        score_matrix[(i,0)] = d * i          # adds i * penalty in Column 0
    
    for j in range(len(string2)+1):
        score_matrix[(0,j)] = d * j          # adds j * penalty in Row 0
        
    for a in range(1, len(string1)+1):
        for b in range(1, len(string2)+1):
            match = score_matrix[(a-1, b-1)] + score[string1[a-1]][string2[b-1]]
            delete = score_matrix[(a-1,b)] + d
            insert = score_matrix[(a,b-1)] + d
            score_matrix[(a,b)] = max(match, insert, delete)         # fills score_matrix with the maximum value between match, mismatch and indel score  
            
    alig_a = ""
    alig_b = ""
    i = len(string1)
    j = len(string2)
    final_score = 0
    
    while i > 0 or j > 0:                # checks score_matrix and creates alignments based on that matrix
         if i > 0 and j > 0 and score_matrix[(i,j)] == score_matrix[(i-1,j-1)] + score[string1[i-1]][string2[j-1]]:
             alig_a = string1[i-1] + alig_a
             alig_b = string2[j-1] + alig_b
             final_score += score[string1[i-1]][string2[j-1]]
             i-= 1
             j-= 1
            
         elif i > 0 and score_matrix[(i,j)] == score_matrix[(i-1,j)] + d:
             alig_a = string1[i-1] + alig_a
             alig_b = "-" + alig_b
             i -= 1
             final_score +=  d
             
         else:
             alig_a = "-" + alig_a
             alig_b = string2[j-1] + alig_b
             j -= 1
             final_score +=  d
     
    return final_score, alig_a, alig_b
```

**Example**\
*Input:*\
PRTEIN\
AKIN\
-1\
BLOSUM62 scoring matrix

*Output:*\
9, 'PRTEIN', '--AKIN'
