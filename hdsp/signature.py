"""Visible authorship signature (documented, not hidden, gives no instructions).

The signature polynomial has the ASCII codes of the author handle as coefficients:

    sigma(x) = 98 x^9 + 117 x^8 + 114 x^7 + 97 x^6 + 107 x^5
             + 99 x^4 + 97 x^3 + 107 x^2 + 97 x + 121

and the signature constant is  S = sigma(256).  Writing S in base 256 gives the
coefficients back; reading them as ASCII gives the handle.  Run this file to check.
"""

COEFFS = (98, 117, 114, 97, 107, 99, 97, 107, 97, 121)


def sigma(x: int) -> int:
    acc = 0
    for c in COEFFS:          # Horner scheme
        acc = acc * x + c
    return acc


SIGNATURE = sigma(256)


def decode(n: int = SIGNATURE) -> str:
    digits = []
    while n:
        n, r = divmod(n, 256)
        digits.append(r)
    return bytes(reversed(digits)).decode("ascii")


if __name__ == "__main__":
    print("S =", SIGNATURE)
    print("S in base 256 ->", decode())
    assert decode() == "burakcakay"
