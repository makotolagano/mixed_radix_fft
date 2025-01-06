import numpy as np

def mixed_radix_digit_reverse(arr, radices):
    """
    Perform digit-reversal for a mixed-radix FFT.

    Parameters:
        arr (np.ndarray): Input array to reorder.
        radices (list[int]): List of radices defining the mixed-radix structure.

    Returns:
        np.ndarray: Reordered array based on digit-reversal.
    """
    n = arr.size
    if n != np.prod(radices):
        raise ValueError("The size of the array must match the product of the radices.")

    def mixed_radix_index(index, radices):
        """Convert a linear index to its mixed-radix representation."""
        digits = []
        for radix in radices:
            digits.append(index % radix)
            index //= radix
        return digits

    def linear_index_from_mixed_radix(digits, radices):
        """Convert a mixed-radix representation back to a linear index."""
        index = 0
        for digit, radix in zip(reversed(digits), reversed(radices)):
            index = index * radix + digit
        return index

    reordered = np.empty_like(arr)
    for i in range(n):
        # Get the mixed-radix representation
        digits = mixed_radix_index(i, radices)
        # Reverse the mixed-radix digits
        reversed_digits = digits[::-1]
        # Map back to a linear index
        new_index = linear_index_from_mixed_radix(reversed_digits, radices)
        reordered[new_index] = arr[i]

    return reordered

def reverse_mixed_radix_digit_reverse(arr, radices):
    """
    Reverse digit-reversal for a mixed-radix FFT.

    Parameters:
        arr (np.ndarray): Input array that was reordered via digit-reversal.
        radices (list[int]): List of radices defining the mixed-radix structure.

    Returns:
        np.ndarray: Array restored to its original order.
    """
    n = arr.size
    if n != np.prod(radices):
        raise ValueError("The size of the array must match the product of the radices.")

    def mixed_radix_index(index, radices):
        """Convert a linear index to its mixed-radix representation."""
        digits = []
        for radix in radices:
            digits.append(index % radix)
            index //= radix
        return digits

    def linear_index_from_mixed_radix(digits, radices):
        """Convert a mixed-radix representation back to a linear index."""
        index = 0
        for digit, radix in zip(reversed(digits), reversed(radices)):
            index = index * radix + digit
        return index

    restored = np.empty_like(arr)
    for i in range(n):
        # Get the mixed-radix representation (from reversed radices)
        digits = mixed_radix_index(i, radices[::-1])
        # Reverse the digits to restore the original index
        original_digits = digits[::-1]
        # Map back to a linear index
        original_index = linear_index_from_mixed_radix(original_digits, radices)
        restored[original_index] = arr[i]

    return restored


data = np.array([0, 1, 2, 3, 4, 5])
radices = [3, 2]
digit_reversed_data = mixed_radix_digit_reverse(data, radices)
print("Digit-Reversed Data:", digit_reversed_data)
restored_data = reverse_mixed_radix_digit_reverse(digit_reversed_data, radices)
print("Restored Original Data:", restored_data)
