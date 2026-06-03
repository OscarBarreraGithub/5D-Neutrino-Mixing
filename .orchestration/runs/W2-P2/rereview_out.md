Independent first-6 roots: 2.450509663813728, 5.567547286728074, 8.701957144582060, 11.840260538783498, 14.980018189299921, 18.120466878759636  
Independent a, bracketed 1024-root tower, Q=2048, N=32 criterion: a(0.2)=21.860547338022837, a(0.65)=-1.457139933034891  
1. RESOLVED: roots are strict, unique, ordered, and near pi-spaced; old n6 skip is gone.  
2. RESOLVED: corrected-tower a(c) matches expected values.  
3. RESOLVED: tests assert strict/unique tower and recompute overlap with independent bracketing.  
4. RESOLVED: endpoint identities, universal subtraction, b_R sign, cache determinism, truncation criterion remain covered/passing.  
pytest `tests/ -q`: 1645 passed, 1 skipped in 807.66s.  
PHASE2-OK