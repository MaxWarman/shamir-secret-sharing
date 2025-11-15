#!/usr/bin/python3

import secrets
from Crypto.Util.number import bytes_to_long, long_to_bytes

from mod_poly import PolyGF

class SSS:
    def __init__(self, total_shares, min_shares, secret, prime):
        self.total_shares = total_shares
        self.min_shares = min_shares
        self.secret = secret
        self.p = prime
        self.poly = self.new_polynomial()

    def new_polynomial(self):
        new_p = self.p
        new_coefficients = [bytes_to_long(self.secret)]
        for _ in range(self.min_shares - 2):
            new_coefficients.append(secrets.randbelow(new_p))
        # Make sure largest exponent coefficient is non-zero
        new_coefficients.append(secrets.randbelow(new_p - 1) + 1)
        return PolyGF(new_coefficients, new_p)
    
    def gen_shares(self):
        shares = []
        for i in range(self.total_shares):
            share_x = i + 1
            share_y = self.poly.eval(share_x)
            shares.append((share_x, share_y))
        return shares

    @classmethod
    def recover_secret(cls, shares, prime):
        recovered_poly = SSS.mod_lagrange_interpolation(shares, prime)
        secret = long_to_bytes(recovered_poly.coefficients[0])
        return secret

    @classmethod
    def mod_lagrange_interpolation(cls, points, prime):
        result_poly = PolyGF.zero(prime)

        for i in range(len(points)):
            current_poly = PolyGF.one(prime)
            for j in range(len(points)):
                if i == j: continue
                current_poly = current_poly * PolyGF([-points[j][0], 1], prime)
            
            denominator = 1
            for j in range(len(points)):
                if i == j: continue
                denominator = (denominator * (points[i][0] - points[j][0])) % prime
            
            current_poly = current_poly * points[i][1] * pow(denominator, -1, prime)
            result_poly = result_poly + current_poly

        return result_poly

def main():
    
    nbit = 256
    SECRET = b"my super secret message!!"
    assert len(SECRET) * 8 < nbit
    
    from Crypto.Util.number import getPrime
    PRIME = getPrime(nbit)
    MIN_SHARES = 3
    TOTAL_SHARES = 7
    sss = SSS(total_shares=TOTAL_SHARES, min_shares=MIN_SHARES, secret=SECRET, prime=PRIME)
    
    assert SECRET == long_to_bytes(sss.poly.coefficients[0])
    
    shares = sss.gen_shares()
    
    from random import sample
    tests = 10
    for _ in range(tests):
        test_shares = sample(shares, MIN_SHARES)
        recovered = SSS.recover_secret(test_shares, PRIME)
        assert recovered == SECRET


if __name__ == "__main__":
    main()