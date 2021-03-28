# Course 3 in the Bioinformatics Specialization: Comparing Genes, Proteins and Genomes

**Week 2 - Local Alignment**

This code is based on the Smith-Waterman algorithm. Usign a ```scoring matrix``` for matches/mismatches and a ```penalty``` variable for insertions or deletions, the function returns the highest possible score for a **local** alignment between the two strings along with what this alignment is.

```python
def local_alignment(string1, string2, penalty, score):     
    d = penalty
    score_matrix = {}
    highest_score = 0
    location = 0
    backtrack = {}
    
    # this function assumes that, up to the point where the alignment begins and after it ends, no score is calculated.
    # for the backtrack matrix "s" will be used to indicate the point in which the alignment ends, "b" is for insertion, "r" for deletion and "c" for matchces/mismatches
    
    for i in range(len(string1)+1):
        score_matrix[(i,0)] = 0          # adds 0 to Column 0
        backtrack[(i,0)] = "s"           # adds "s" to Column 0 
    
    for j in range(len(string2)+1):
        score_matrix[(0,j)] = 0          # adds 0 to Row 0
        backtrack[(0,j)] = "s"           # adds "s" to Row 0
        
    for a in range(1, len(string1)+1):
        for b in range(1, len(string2)+1):
            match = score_matrix[(a-1, b-1)] + score[string1[a-1]][string2[b-1]]
            delete = score_matrix[(a-1,b)] + d
            insert = score_matrix[(a,b-1)] + d
            score_matrix[(a,b)] = max(0, match, insert, delete)      # based on the scoring matrix provided, fills the score_matrix with maximum value between match, mismatch and indel penalty
            
            if score_matrix[(a,b)] > highest_score:                  # keeps the highest score in score_matrix, this will be the starting point for backtracking
                highest_score = score_matrix[(a,b)]
                location = (a,b)
            
            if score_matrix[(a,b)] == 0:                             # fills backtrack matrix based on score_matrix
                backtrack[(a,b)] = "s"
            elif score_matrix[(a,b)] == delete:
                backtrack[(a,b)] = "d"
            elif score_matrix[(a,b)] == insert:
                backtrack[(a,b)] = "r"
            else:
                backtrack[(a,b)] = "c"
                
    alig_a = ""
    alig_b = ""
    i = location[0]
    j = location[1]
    final_score = 0
    
    while i > 0 or j > 0:                # checks backtrack matrix and creates alignment starting at the point with max score stopping when it reaches an "s" point
        
        if backtrack[(i,j)] == "c":
            alig_a = string1[i-1] + alig_a
            alig_b = string2[j-1] + alig_b
            i -= 1
            j -= 1
            
        elif backtrack[(i,j)] == "r":
            alig_a = "-" + alig_a
            alig_b = string2[j-1] + alig_b
            j -= 1
            
        elif backtrack[(i,j)] == "d":
            alig_a = string1[i-1] + alig_a
            alig_b = "-" + alig_b
            i -= 1
            
        else:
            break
    
    for letter in range(len(alig_a)):     # calculates de score of the local alignment created to check if it is the same as the max score
        if alig_a[letter] == "-":
            final_score += d
        elif alig_b[letter] == "-":
            final_score += d
        else:
            final_score += score[alig_a[letter]][alig_b[letter]]
                
    return final_score, highest_score, alig_a, alig_b
```

**Example**\
*Input:*\
PRTEIN\
AKIN\
-1\
BLOSUM62 scoring matrix

*Output:*\
11, 11, 'EIN', 'KIN'
