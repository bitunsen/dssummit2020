import unittest
import random

from EightQueenProblem import Gene, Chromosome, GAContext, generate_initial_population, \
    execute_generate_eight_queen_solution


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

    def test_chromosome_crossover(self):
        _data_array_1 = [i for i in range(8)]
        random.shuffle(_data_array_1)
        _parent_1 = Chromosome(_data_array_1)

        _data_array_2 = [i for i in range(8)]
        random.shuffle(_data_array_2)
        _parent_2 = Chromosome(_data_array_2)

        try:
            _child_1, _child_2 = _parent_1.crossover(_parent_2)

            if _child_1 == _parent_1 or _child_1 == _parent_2:
                self.fail("Test case failed - Crossover is producing children with same as one of the parents")
            elif _child_2 == _parent_1 or _child_2 == _parent_2:
                self.fail("Test case failed - Crossover is producing children with same as one of the parents")

        except Exception:
            self.fail("Unexpected Error occurred")
        else:
            pass

    def test_chromosome_mutation(self):
        _data_array_1 = [i for i in range(8)]
        random.shuffle(_data_array_1)
        _parent_1 = Chromosome(_data_array_1)

        try:
            _child_1 = _parent_1.mutation(0.900)

            if _child_1 == _parent_1:
                self.fail("Test case failed - Mutation is producing children with same as the parent")
            else:
                pass

        except Exception:
            self.fail("Unexpected Error occurred")
        else:
            pass

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

    def test_final_result(self):
        try:
            _target_chromosome, is_target_solution_found = execute_generate_eight_queen_solution(8, 10, 0.1)

            if is_target_solution_found:
                print("The target sequence {}.".format(_target_chromosome.get_raw_data_sequence()))
                pass
            else:
                self.fail("Unable to find any target sequence withing 2000 generations.")

        except (TypeError, ValueError):
            self.fail("Unexpected error occurred")


if __name__ == '__main__':
    unittest.main()
