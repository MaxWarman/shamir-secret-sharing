#!/usr/bin/python3

from Crypto.Util.number import isPrime

class PolyGF:
    def __init__(self, coefficients, prime):
        if not len(coefficients) > 0: raise ValueError("PolyGF: Polynomial must have at least one coefficient!")
        if not isPrime(prime): raise ValueError(f"PolyGF: Provided number is not prime: {prime}")
        # if coefficients[-1] % prime == 0: raise ValueError(f"PolyGF: Highest exponent's coefficient must be non-zero: {coefficients[-1]} % {prime} = 0")
        # for c in coefficients:
        #     if not 0 <= c < prime: raise ValueError(f"PolyGF: Coefficients must be in range [0, {prime}): {c}")
        
        for i in range(len(coefficients)):
            coefficients[i] %= prime

        self.coefficients = list(coefficients)
        self.p = prime
        self.reduce()
    
    def reduce(self):
        inv_coefficients = list(self.coefficients)[::-1]
        while len(inv_coefficients) > 1:
            if inv_coefficients[0] % self.p != 0:
                break
            inv_coefficients.pop(0)
        self.coefficients = list(inv_coefficients)[::-1]
        self.deg = (len(self.coefficients) - 1)

    @classmethod
    def add(cls, poly1, poly2):
        if poly1.p != poly2.p: raise ValueError(f"PolyGF.add: Added polynomials must have the same field order: {poly1.p}, {poly2.p}")
        
        result_p = poly1.p
        result_deg = max(len(poly1.coefficients), len(poly2.coefficients))

        result_coefficients = [0 for _ in range(result_deg + 1)]

        for i in range(len(result_coefficients)):
            if i < len(poly1.coefficients):
                result_coefficients[i] += poly1.coefficients[i]
            if i < len(poly2.coefficients):
                result_coefficients[i] += poly2.coefficients[i]
            result_coefficients[i] %= result_p

        return PolyGF(result_coefficients, result_p)
    
    def __add__(self, other):
        if not isinstance(other, PolyGF): raise ValueError(f"PolyGF.__add__: Added elements must be polynomials: {type(other)}")
        return PolyGF.add(self, other)
    
    @classmethod
    def sub(cls, poly1, poly2):
        if poly1.p != poly2.p: raise ValueError(f"PolyGF.sub: Subtracted polynomials must have the same field order: {poly1.p}, {poly2.p}")
        
        result_p = poly1.p
        result_deg = max(len(poly1.coefficients), len(poly2.coefficients))
        result_coefficients = [0 for _ in range(result_deg + 1)]

        for i in range(len(result_coefficients)):
            if i < len(poly1.coefficients):
                result_coefficients[i] += poly1.coefficients[i]
            if i < len(poly2.coefficients):
                result_coefficients[i] -= poly2.coefficients[i]
            result_coefficients[i] %= result_p

        return PolyGF(result_coefficients, result_p)
    
    def __sub__(self, other):
        if not isinstance(other, PolyGF): raise ValueError(f"PolyGF.__sub__: Subtracted elements must be polynomials: {type(other)}")
        return PolyGF.sub(self, other)
    
    def eval(self, argument):
        value = 0
        for i, c in enumerate(self.coefficients):
            value = (value + c * pow(argument, i, self.p)) % self.p
        return value

    @classmethod
    def mul_poly(cls, poly1, poly2):
        if poly1.p != poly2.p: raise ValueError(f"PolyGF.mul: Multiplied polynomials must have the same field order: {poly1.p}, {poly2.p}")
        
        result_p = poly1.p
        if PolyGF.eq(poly1, PolyGF.zero(result_p)) or PolyGF.eq(poly2, PolyGF.zero(result_p)): return PolyGF.zero(result_p)
        if PolyGF.eq(poly1, PolyGF.one(result_p)): return PolyGF(poly2.coefficients, poly2.p)
        if PolyGF.eq(poly2, PolyGF.one(result_p)): return PolyGF(poly1.coefficients, poly1.p)
        
        result_coefficients = [0 for _ in range(poly1.deg + poly2.deg + 1)]

        for i, c1 in enumerate(poly1.coefficients):
            for j, c2 in enumerate(poly2.coefficients):
                result_coefficients[i + j] = (result_coefficients[i + j] + (c1 * c2)) % result_p

        return PolyGF(result_coefficients, result_p)
    
    @classmethod
    def mul_scalar(cls, poly1, scalar):
        result_coefficients = []
        result_p = poly1.p
        for coef in poly1.coefficients:
            new_coefficient = (coef * scalar) % result_p
            result_coefficients.append(new_coefficient)

        return PolyGF(result_coefficients, result_p)
    
    def __mul__(self, other):
        if isinstance(other, PolyGF): return PolyGF.mul_poly(self, other)
        elif isinstance(other, int): return PolyGF.mul_scalar(self, other)
        else: raise ValueError(f"PolyGF.__mul__: Cannot handle multiplication by {type(other)} type.")

    def __rmul__(self, other):
        if isinstance(other, PolyGF): return PolyGF.mul_poly(self, other)
        elif isinstance(other, int): return PolyGF.mul_scalar(self, other)
        else: raise ValueError(f"PolyGF.__mul__: Cannot handle multiplication by {type(other)} type.")

    @classmethod
    def eq(cls, poly1, poly2):
        if poly1.p != poly2.p: raise ValueError(f"PolyGF.eq: Compared polynomials must have the same field order: {poly1.p}, {poly2.p}")
        poly1.reduce()
        poly2.reduce()
        if list(poly1.coefficients) == list(poly2.coefficients):
            return True
        return False
    
    @classmethod
    def zero(cls, prime):
        return PolyGF([0], prime)
    
    @classmethod
    def one(cls, prime):
        return PolyGF([1], prime)
    
    def __str__(self):
        txt = f"{self.coefficients}; deg: {self.deg}; p: {self.p}"
        return txt


def main():
    print("main")

if __name__ == "__main__":
    main()