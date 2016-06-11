"################################# Classes ##################################"


class Coder(object):

    def __init__(
            self,
            length=6,
            alphabet=(
                "abcdefghijklmnopqrstuvwxyz"
                "123456789"
                "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                ),
            ):
        self.length = length
        self.alphabet = alphabet
        self.n = len(self.alphabet)
        return

    def __getitem__(self, key):
        return self.index(key)

    def index(self, pattern):
        if len(pattern) == 0:
            return 0

        if isinstance(pattern, str):
            pattern = list(pattern)

        character = pattern.pop(-1)
        return self.n * self.index(pattern) + self.alphabet.index(character)

    def pattern(self, index, k=None):
        if k is None:
            k = self.length

        assert k > 0
        if k == 1:
            return self.alphabet[index]

        prefix_index = index // self.n
        r = index % self.n
        return self.pattern(prefix_index, k - 1) + self.alphabet[r]
