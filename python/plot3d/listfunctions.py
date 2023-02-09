
def unique_pairs(listOfItems=[]):
    """Checks if an item is not already in the list 
    
    Args:
        listOfItems (list): list of combinations e.g. (1,2),(3,4),(2,1)

    Yields:
        unique pair: [description]
    """
    seen = set()  #use set to keep track of already seen items, sets provide O(1) lookup  
    for x,y in listOfItems:
        if x!=y and (y,x) not in seen:
            seen.add((x,y)) 
            yield x,y
