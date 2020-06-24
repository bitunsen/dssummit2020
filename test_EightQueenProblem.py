import unittest

from EightQueenProblem import Gene, Chromosome, GAContext, generate_initial_population


def is_valid_chromosome(chromosome, num_of_queens):
    is_valid = True
    for _data in chromosome.get_raw_data_sequence():
        if not isinstance(_data, int):
            is_valid = False
        if 0 > _data > num_of_queens:
            is_valid = False
    return is_valid

class EightQueenProblemTestCase(unittest.TestCase):
    def test_gene_constructor(self):
        try:
            gene_1 = Gene(4, "5")
        except TypeError:
            pass
        else:
            self.fail("Unknown Exception raised.")

    def test_gene_equals_true(self):
        gene_1 = Gene(4, 5)
        gene_2 = Gene(4, 5)
        expected = gene_1.equals(gene_2)
        self.assertEqual(expected, True)

    def test_gene_equals_false(self):
        gene_1 = Gene(4, 5)
        gene_2 = Gene(4, 8)
        expected = gene_1.equals(gene_2)
        self.assertEqual(expected, False)

    def test_gene_equals_type_check(self):
        gene_1 = Gene(4, 5)
        self.assertRaises(TypeError, gene_1.equals, "String")

    def test_gene_prints_type_check(self):
        gene_1 = Gene(4, 5)
        self.assertRaises(TypeError, gene_1.print, 123)

    def test_gene_prints_type_check(self):
        gene_1 = Gene(4, 5)
        self.assertRaises(TypeError, gene_1.print, 123)

    def test_chromosome_constructor_invalid_element_type(self):
        try:
            Chromosome(['a', 'b', 'c', 'd'])
        except TypeError:
            pass
        else:
            self.fail("Unknown Exception raised.")

    def test_chromosome_constructor_negative_integer(self):
        try:
            Chromosome([-1, -2, 4, 5, 6])
        except TypeError:
            pass
        else:
            self.fail("Unknown Exception raised.")

    def test_chromosome_constructor_valid_data(self):
        try:
            Chromosome([0, 1, 2, 3, 4, 5, 6])
            pass
        except TypeError:
            self.fail("TypeError exception raised with valid data.")

    def test_chromosome_constructor_valid_data_invalid_length(self):
        try:
            Chromosome([0])
        except ValueError:
            pass
        else:
            self.fail("Unknown Exception raised.")



    def test_generate_initial_population(self):
        number_of_queens = 8
        execution_context = GAContext(10, number_of_queens, 0.1)
        _initial_population_pool = generate_initial_population(execution_context)
        _is_valid = True
        for _chromosome in _initial_population_pool:
            _is_valid = is_valid_chromosome(_chromosome, number_of_queens)
            if not _is_valid:
                break
        if not _is_valid:
            raise ValueError("Test Failed - ")


if __name__ == '__main__':
    unittest.main()
