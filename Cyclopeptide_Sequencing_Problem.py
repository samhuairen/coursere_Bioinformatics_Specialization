def compute_subpepetides_number(string_length):
    """ return the number of all subpepetides given the length of the pepetides"""
    if isinstance(string_length, int):
        return string_length * (string_length - 1)
    else:
        print("you should input a integer")
    return

