"""
Created 16 Oct 2016.

Description: 
A simple script that returns the list of perfect numbers which are inferior to a given number n.

@author: ghaskan
"""

# This Python file uses the following encoding: utf-8


def divisors(n): # auxiliary function
    """
    :param n: Integer for which we want to find divisors of.
    :return: The list of divisors of n, except for itself.
    """
    if type(n) != int:
        raise TypeError("Please input an integer.")
    d = []
    for i in range(1,n):
        if n % i == 0:
            d.append(i)
    return d


def perfect_numbers(n): # main function
    """
    :param n: An integer.
    :return: The list of perfect numbers inferior to n.
    """
    if type(n) != int:
        raise TypeError("Please input an integer.")
    if n < 0:
        raise ValueError("n cannot be negative.")
    if n <= 6:
        return '{} {}'.format("There are no perfect numbers below", n)
    perfect = []
    for i in range(1, n):
        d = divisors(i)
        if i == sum(d):
            perfect.append(i)
    return '{} {} {} {}'.format("The list of perfect numbers below", n,
                                "are:", ', '.join([str(x) for x in perfect]))
