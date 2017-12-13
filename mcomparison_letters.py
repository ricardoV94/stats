import numpy as np
from string import ascii_lowercase
from itertools import combinations

def mcomparison_letters(p_values, significance_level = 0.05, return_numbers = False):
    """
    Compute letters that represent multiple comparisons between groups
   
    Arguments:
    ---------------
    p_valies: float array
        corrected p-values for each pairwise comparison,
        e.g. for 3 groups, [0.03, 0.001, 0.6] to indicate that at 0.05
        significance level: group1 != group 2; group1 != group3; and group2 = group3

    significance_level: float
        error rate used to establish significancy of differences

    return_numbers: bool
        if True array of numbers is returned, instead of letters

    Returns:
    ---------------
    final_letters: list of lists
        list of list of letters/ numbers that describe multiple comparisons for each group
        e.g. from above, [['a'], ['b'], ['b']]
        e.g. from above for significance_level of 0.01, [['a'], ['a', 'b'], ['b']]
    """
        
    
    # Input preparation
    # -------------------------------------------------------------------------
    
    p_values = np.array(p_values).flatten()    
    
    # Calculate number of groups
    n_groups = 0
    comparisons = len(p_values)
    while comparisons > 0:
        n_groups += 1
        comparisons -= n_groups
    n_groups += 1
    if comparisons < 0:
        raise ValueError('Invalid number of p_values. Must be a triangual number (1+2+3+4+...+n)')
            
    # Convert p-values to boolean and 
    # Rearrange them into a nested list [[1vs2, 1vs3, 1vs4], [2vs3, 2vs4], [3vs4]]
    p_values = p_values <= significance_level
    nested_p_values = [[] for _ in range(n_groups - 1)]
    for icomb, comb in enumerate((combinations(range(n_groups), 2))):
        nested_p_values[comb[0]].append(p_values[icomb])
            
        
    # Letter assignment loop
    # ------------------------------------------------------------------------- 
    
    final_letters = [set([]) for _ in range(n_groups)]
    current_letter, max_letter, new_value_assigned = 1, 1, False
    for i in range(n_groups):                        
        new_value_assigned = False
                
        # If set is empty: Automatically assign new letter
        if not len(final_letters[i]):
            final_letters[i].add(current_letter)
            new_value_assigned = True
            
        if i == len(final_letters)-1:
            break
                 
        # Loop through next groups, looking for similar cases
        current_letters = final_letters[i]
        for j, significant_pval in enumerate(nested_p_values[i]):
            
            if not significant_pval:
            
                # Check if similar group already shares the letter
                for letter in final_letters[i+j+1]:
                    if letter in current_letters:
                        break
                        
                # If it does not, must assign it
                else:
                    new_value_assigned = True
                    offset = 0                   
                    
                    # Check if there are conflicting groups in between that may require adding a new letter
                    # Example instead of (i: a, mid: a, j: a), we have (i: ab, mid: a, j:b)
                    for mid, significant_mid_pval in enumerate(nested_p_values[i][:j]):
                        if not significant_mid_pval: #(i = mid)             
                            if nested_p_values[i+mid+1][j-mid-1]: #(mid != j)    
                                if current_letter+offset in final_letters[i+mid+1]: #(mid contains new letter)
                                # In the example: [['a', 'b'], ['a', 'c'], ['a'], ['b', 'c']]
                                # The 'b' (1st loop) is a result of an offset, whereas the 'c' (2nd loop) is not
                                # because the 3rd group is already assigned a different value (i.e., 'a') by then  
                                    offset += 1
                                                                        
                    # If there is a new letter assing it retroactively (i:a -> i:a,b)
                    # ... And finally assing it to the new group (j: new letter)
                    final_letters[i].add(current_letter+offset)
                    final_letters[i+j+1].add(current_letter+offset)
                        
                    if current_letter+offset > max_letter:
                        max_letter = current_letter+offset  
                        
        # Update to next letter
        if new_value_assigned and max_letter >= current_letter:
            max_letter += 1
            current_letter = max_letter


    # Return output
    # -------------------------------------------------------------------------
    
    if return_numbers:
        return [[l for l in letter ] for letter in final_letters]
    
    alphabet = dict(enumerate(ascii_lowercase, 1))
    return [[alphabet[l] for l in letter ] for letter in final_letters]