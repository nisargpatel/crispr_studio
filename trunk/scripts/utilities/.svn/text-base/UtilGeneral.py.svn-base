#!/usr/bin/env python

"""
This module provides various computational utility functions.
"""

__author__ = "Ruben Acuna"
__copyright_ = "Copyright(c) 2011, ASU iGEM Team"

def SAT(expression, text):
    """Checks if the given text is satisfied by the given boolean expression. The grammar we use is listed within."""

    if expression == "":
        return True

    expression = expression.lower()
    text = text.lower()

    terms = expression.split(" or ")
    matches = False

    #We should support the following:
    #   <expression> ::= <term> (OR <term>)*
    #   <term>       ::= <factor> [AND <factor>]*
    #   <factor>     ::= [NOT] <literal> | <literal>

    #HACK: need to use LEX and YACC.

    for term in terms:
        termTokens = term.split(" ")
        currentFactor = 0
        termMatches = True

        while currentFactor < len(termTokens):
            if termTokens[currentFactor] == "not":
                #consume the NOT
                currentFactor += 1

                if currentFactor < len(termTokens):
                    if termTokens[currentFactor] in text:
                        termMatches = False

                #consume the literal
                currentFactor += 1

            else:
                if not termTokens[currentFactor] in text:
                    termMatches = False

                #consume the literal
                currentFactor += 1

            #consume the tailing AND
            currentFactor += 1

        if termMatches:
            matches = True

    return matches