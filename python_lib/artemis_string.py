''' Function to handle artemis strings '''

def filter_string(s):
    ''' From Artemis to normal string '''
    return s.replace("%2C", ',').replace("+", " ").replace("%28", "(").replace("%09", "").\
             replace("%29", ")").replace("%3D", "=").replace("%2B", "+").replace("%2F","/").\
             replace("%27", "'").replace("%3B", ";")

def artemis_string(s):
    ''' From normal string to Artemis string '''
    return s.replace(',', "%2C").replace("+", "%2B").replace(" ", "+").\
             replace("(", "%28").replace(")", "%29").replace("=", "%3D")