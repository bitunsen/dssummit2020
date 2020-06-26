import numpy as np
import random

system_random = random.SystemRandom()


def get_random_element_for_mutation(valid_gene_data):
    return system_random.choice(valid_gene_data)


def get_random_number_from_list(number_list):
    return system_random.choice(number_list)


"""
Gene is a class holding the row and column values of a particular position in Chess board
x_cor denotes the row number and y_cor denotes the column number
"""


class Gene:

    def __init__(self, x_cor, y_cor):
        if not isinstance(x_cor, int) or not isinstance(y_cor, int):
            raise TypeError("Gene can take only positive integer values")
        if x_cor < 0 or y_cor < 0:
            raise TypeError("Gene can take only positive integer values")

        self.x_cor = x_cor
        self.y_cor = y_cor

    def equals(self, point):
        if not isinstance(point, Gene):
            raise TypeError("Objects of type Gene is allowed only")

        is_same = False
        if self.x_cor == point.x_cor and self.y_cor == point.y_cor:
            is_same = True
        return is_same

    def print(self, message):
        if not isinstance(message, str):
            raise TypeError("Only String message is allowed.")

        print(message, " (", self.x_cor, ", ", self.y_cor, ")")


"""
Chromosome is a class containing an array of instances of Gene class
An instance of Chromosome will hold the information of the positions of a set of 8 queens
in a chess board. It also holds the fitness score of that solution.

Mutation rate is calculated based on the uniqueness of data elements in the data array
It is calculated by dividing the number of unique elements by the length of data elements

"""


class Chromosome:

    def __init__(self, _input_data_list):
        if not isinstance(_input_data_list, list):
            raise TypeError("Invalid argument passed - only array of positive integers are allowed.")

        if len(_input_data_list) < 2:
            raise ValueError("Invalid argument passed - number of elements must be more than 2 in the input array.")

        for _element in _input_data_list:
            if not isinstance(_element, int):
                raise TypeError("Invalid argument passed - only array of positive integers are allowed.")
            if _element < 0:
                raise TypeError("Invalid argument passed - only array of positive integers are allowed.")

        self._data_list = []
        self._gene_list = []
        self.fitness_score = 0
        self._mutation_score = 0.0

        self._data_list = _input_data_list[:]
        for i in range(len(_input_data_list)):
            _point = Gene(i, _input_data_list[i])
            self._gene_list.append(_point)
        self._mutation_score = random.uniform(0, 1)

    def get_mutation_score(self):
        return self._mutation_score

    def print_chromosome(self, message):
        gene_size = len(self._data_list)
        chromosome_str = "["
        for i in range(gene_size):
            if i < (gene_size - 1):
                chromosome_str = chromosome_str + str(self._data_list[i]) + ", "
            if i == (gene_size - 1):
                chromosome_str = chromosome_str + str(self._data_list[i])
        chromosome_str = chromosome_str + "]"
        print(message, "Chromosome : ", chromosome_str, " ::: Score : ", self.fitness_score,
              " :: Mutation Score :: ", self.get_mutation_score())

    def get_fitness_score(self):
        return self.fitness_score

    def set_fitness_score(self, score):
        self.fitness_score = score

    def get_raw_data_sequence(self):
        return self._data_list

    def get_gene_sequence(self):
        return self._gene_list

    def equals(self, _other_chromosome):
        if not isinstance(_other_chromosome, Chromosome):
            raise TypeError("Invalid argument passed - only object of type Chromosome is allowed.")
        is_equal = True
        for index in range(len(_other_chromosome.get_gene_sequence())):
            _other_gene = _other_chromosome.get_gene_sequence()[index]
            _self_gene = self.get_gene_sequence()[index]
            if not _self_gene.equals(_other_gene):
                is_equal = False
                break
        return is_equal

    """
    Rather than having fixed point crossover for all the parent chromosome, select 
    the crossover point randomly. This introduces more randomness in the chromosome during crossover 
    """

    def crossover(self, other_parent):
        if not isinstance(other_parent, Chromosome):
            raise TypeError("Invalid argument passed - only object of type Chromosome is allowed.")

        if len(other_parent.get_raw_data_sequence()) != len(self.get_raw_data_sequence()):
            raise TypeError(
                "Invalid argument passed - chromosome passed as argument does not have same number of Gene.")

        index_crossover_array = [i for i in range(len(self.get_raw_data_sequence()))]
        index_crossover = get_random_number_from_list(index_crossover_array)
        if index_crossover == 0 or index_crossover == (len(self.get_raw_data_sequence()) - 1):
            index_crossover = int(len(self.get_raw_data_sequence()) / 2)

        # self.print_chromosome("CROSSOVER :: Parent 1 Chromosome :: ")
        # other_parent.print_chromosome("CROSSOVER :: Parent 1 Chromosome :: ")
        # print("CROSSOVER ::: Selected Crossover Index :: ", index_crossover)

        parent_1_data = self.get_raw_data_sequence()
        parent_2_data = other_parent.get_raw_data_sequence()
        child_data_1 = parent_1_data[0:index_crossover]
        child_data_2 = parent_2_data[0:index_crossover]

        set_1 = set(child_data_1)
        set_2 = set(child_data_2)

        for i in range(index_crossover, len(self.get_raw_data_sequence())):
            for j in range(len(self.get_raw_data_sequence())):
                if not parent_2_data[j] in set_1:
                    child_data_1.append(parent_2_data[j])
                    set_1.add(parent_2_data[j])

                if not parent_1_data[j] in set_2:
                    child_data_2.append(parent_1_data[j])
                    set_2.add(parent_1_data[j])

        # print("AFTER CROSSOVER :: Child 1 Data :: ", child_data_1)
        # print("AFTER CROSSOVER :: Child 2 Data :: ", child_data_2)
        _random_sample_index_array = random.sample(range(len(child_data_1)), 2)
        _index_1 = _random_sample_index_array[0]
        _index_2 = _random_sample_index_array[1]
        child_data_1[_index_1], child_data_1[_index_2] = child_data_1[_index_2], child_data_1[_index_1]

        _random_sample_index_array = random.sample(range(len(child_data_2)), 2)
        _index_1 = _random_sample_index_array[0]
        _index_2 = _random_sample_index_array[1]
        child_data_2[_index_1], child_data_2[_index_2] = child_data_2[_index_2], child_data_2[_index_1]

        child_1 = Chromosome(child_data_1)
        child_2 = Chromosome(child_data_2)

        """
        Being Random in nature, sometimes parents and children are becoming same, and as a result
        the whole process is going to infinite loop. This below code is the circuit-breaker.
        If parent and child are same, then do a shuffle of the child
        """
        if self.equals(child_1) or other_parent.equals(child_1):
            random.shuffle(child_data_1)
            child_1 = Chromosome(child_data_1)

        if self.equals(child_2) or other_parent.equals(child_2):
            random.shuffle(child_data_2)
            child_2 = Chromosome(child_data_2)

        return child_1, child_2

    """
    This method selects two index numbers randomly. Then it swaps the value of the elements.
    If the input mutation rate is greater than the random number assigned to this chromosome
    then perform the mutation
    """

    def mutation(self, input_mutation_rate):

        if self._mutation_score < input_mutation_rate:
            if len(self._data_list) == 2:
                _index_1 = 1
                _index_2 = 0
            else:
                _random_sample_index_array = random.sample(range(len(self._data_list)), 2)
                _index_1 = _random_sample_index_array[0]
                _index_2 = _random_sample_index_array[1]

            self._data_list[_index_1], self._data_list[_index_2] = self._data_list[_index_2], self._data_list[_index_1]
            point = Gene(_index_1, self._data_list[_index_1])
            self._gene_list[_index_1] = point
            point = Gene(_index_2, self._data_list[_index_2])
            self._gene_list[_index_2] = point


"""
This utility method takes the value of one position of chess board 
and returns the all possible attacking positions 
"""


def get_all_attacking_points_from_one_point(x_cor, y_cor, num_of_queens):
    if not isinstance(x_cor, int) or not isinstance(y_cor, int) or not isinstance(num_of_queens, int):
        raise TypeError("Invalid argument passed - Only positive integer is allowed.")

    if 0 > x_cor < num_of_queens or 0 > y_cor < num_of_queens:
        raise ValueError("Invalid argument passed - Row and Column values must "
                         "be an integer less than number of queens.")

    _invalid_point_arrays = []

    for i in range(num_of_queens):
        if i == x_cor:
            continue
        offset = (x_cor - i)

        if 0 <= (y_cor + offset) < num_of_queens:
            _invalid_point_arrays.append(Gene(i, (y_cor + offset)))

        if 0 <= (y_cor - offset) < num_of_queens:
            _invalid_point_arrays.append(Gene(i, (y_cor - offset)))

    for i in range(num_of_queens):
        if i == x_cor:
            continue
        _invalid_point_arrays.append(Gene(i, y_cor))
        # print("Point :: (",x_cor,", ",i,")")
    return _invalid_point_arrays


class GAContext:

    def __init__(self, size_of_population, num_of_queens, mutation_rate):
        if not isinstance(size_of_population, int):
            raise TypeError("Invalid input provided - Size of population must of a positive Integer")

        if not isinstance(num_of_queens, int):
            raise TypeError("Invalid input provided - Number of Queens must of a positive Integer")

        if num_of_queens != 8:
            raise ValueError("WARNING - This program is not optimized to handle number of queens greater than 8.")

        self._attacking_point_dictionary = {}
        self._valid_gene_data = []
        self._population_pool = []
        self._max_generation = 3000
        self._system_random = random.SystemRandom()
        self._current_generation = 0
        self._size_of_population = size_of_population
        self._mutation_rate = mutation_rate

        self._length_of_chromosome = num_of_queens
        self._valid_gene_data = [i for i in range(num_of_queens)]
        self.__populate_attacking_point_dictionary()

    def get_mutation_rate(self):
        return self._mutation_rate

    def get_size_of_population(self):
        return self._size_of_population

    def get_valid_gene_data(self):
        return self._valid_gene_data

    def get_current_generation(self):
        return self._current_generation

    def set_current_generation(self, cur_generation):
        self._current_generation = cur_generation

    def get_num_of_queens(self):
        return self._length_of_chromosome

    def __populate_attacking_point_dictionary(self):
        num_of_queens = self._length_of_chromosome
        for row in range(num_of_queens):
            for column in range(num_of_queens):
                key = str(row) + "-" + str(column)
                self._attacking_point_dictionary[key] = get_all_attacking_points_from_one_point(row, column,
                                                                                                num_of_queens)

    def get_all_attacking_points_from_dictionary(self, gene_element):
        x_cor, y_cor = gene_element.x_cor, gene_element.y_cor
        key = str(x_cor) + "-" + str(y_cor)
        return self._attacking_point_dictionary[key]

    def get_population_pool_size(self):
        return len(self._population_pool)

    def get_population_pool(self):
        return self._population_pool

    def initialize_population_pool(self, _input_population_pool):
        self._population_pool = []
        self._population_pool = _input_population_pool[:]

    # Select parents using tournament process
    # Using Tournament process, Only half of the parents will be selected for reproduction
    # This method will not be used as tournament based parent selection didn't give good results
    def run_parent_selection_process_tournament(self):
        _max_score = 0
        _score_array = []
        _selected_parent_pool = []
        _selected_parent_index_set = set()
        _selected_parent_pool_size = int(self._size_of_population / 2)

        # print("-------------- START OF TOURNAMENT METHOD OF PARENT SELECTION --------------------------")
        # print("[GEN - ", self._current_generation, "] :: Number of parents to be selected :: ",
        #      _selected_parent_pool_size)

        # Calculate Fitness Score (Number of non-attacking Queens) for each chromosome
        for _chromosome in self._population_pool:
            # _chromosome_score = self.calculate_fitness_score(_chromosome)
            _score_array.append(_chromosome.get_fitness_score())

        if len(self._population_pool) > self._size_of_population:
            self._population_pool = self._population_pool[:self._size_of_population]

        _max_score = np.max(_score_array)

        _count = 0
        while len(_selected_parent_index_set) < _selected_parent_pool_size:
            _random_sample_index = random.sample(range(len(self._population_pool)), 4)
            _count += 1
            if _count % 10 == 0:
                print("Infinite Loop ....")
            _max_selected_score = -1
            _max_selected_index = -1
            for _index in _random_sample_index:
                if self._population_pool[_index].get_fitness_score() > _max_selected_score:
                    _max_selected_index = _index
            if _max_selected_index in _selected_parent_index_set:
                continue
            else:
                _selected_parent_pool.append(self._population_pool[_max_selected_index])
                _selected_parent_index_set.add(_max_selected_index)

        # for _chromosome in _selected_parent_pool:
        #    _chromosome.print_chromosome("[GEN - " + str(self._current_generation) + "] :: Selected Parent :: ")

        # selected_parent_pool = self._population_pool[:]

        return _selected_parent_pool, _max_score


"""
This method generates '_size_of_population' number of Chromosomes
as part of initial population. 
"""


def generate_initial_population(_context):
    _context.set_current_generation(1)
    for count_population in range(_context.get_size_of_population()):
        # print("Valid Data Array :: ", self._valid_gene_data)
        gene_array = _context.get_valid_gene_data()[:]
        random.shuffle(gene_array)
        # print("Gene Array :: ", gene_array)
        _initial_chromosome = Chromosome(gene_array)
        _chromosome_score = calculate_fitness_score(_context, _initial_chromosome)
        _context.get_population_pool().append(_initial_chromosome)

    # Sorting (in ascending order) all the chromosomes present in the population pool
    _context.initialize_population_pool(sorted(_context.get_population_pool(),
                                               key=lambda element: element.fitness_score, reverse=True))
    # for _element in self._population_pool:
    #    _element.print_chromosome("Initial Pool :: ")

    # print("-------------- END OF POPULATION GENERATION AND SORTING BASED ON SCORE --------------------------")
    return _context.get_population_pool()


"""
This method reads a line containing all valid solutions of 8 Queens. It reads each line and 
compare with the input solution

"""


def validate_fitness_score(_data_array):
    _is_valid = False
    with open('8queen_all_solutions.txt', 'r') as file_1:
        for line in file_1:
            line = line.rstrip()
            _str_array = line.split(",")
            _num_array = [int(str_value) for str_value in _str_array]

            if _num_array == _data_array:
                _is_valid = True
                break

    return _is_valid


"""
    This function is responsible to calculate the fitness score of one particular chromosome
    It finds out total number of non-attacking queen pairs for each Gene (position of a queen)
    and then sum it up. Now, one gene can have maximum 7 non-attacking queens, in that case, 
    the fitness score of a correct sequence will be 56 (7*8).  
"""


def calculate_fitness_score(_context, input_chromosome):
    if not isinstance(input_chromosome, Chromosome):
        raise TypeError("Invalid argument - Only object of type Chromosome is expected.")

    _data = input_chromosome.get_raw_data_sequence()[:]
    _non_attacking_queen = [0] * _context.get_num_of_queens()
    for index in range(len(_data)):
        current_point = Gene(index, _data[index])
        # current_point.print("Current Queen :: ")
        attacking_queen = 0
        current_attacking_points = _context.get_all_attacking_points_from_dictionary(current_point)
        for i in range(len(_data)):
            if i == index:
                continue
            other_queen = Gene(i, _data[i])
            for attacking_point in current_attacking_points:
                if other_queen.equals(attacking_point):
                    attacking_queen = attacking_queen + 1
                    # other_queen.print("Attacking Queen :: ")

        _non_attacking_queen[index] = ((len(_data) - 1) - attacking_queen)
        # print(non_attacking_queen)
        input_chromosome.fitness_score = np.sum(_non_attacking_queen)
    return input_chromosome.fitness_score


# This method assumes that the chromosome pool is already sorted
# Take half of the high scored parents for producing next generation kids
def run_parent_selection_process(_context):
    _max_score = 0
    _score_array = []
    _selected_parent_pool = []
    _selected_parent_index_set = set()
    _selected_parent_pool_size = int(_context.get_size_of_population() / 2)

    if _selected_parent_pool_size % 2 != 0:
        _selected_parent_pool_size += 1

    _selected_parent_pool = _context.get_population_pool()[0:_selected_parent_pool_size]
    return _selected_parent_pool, _max_score


"""
This method will perform the crossover and mutation on a selected set of parents
and produce children for next generation
"""


def execute_crossover_mutation(_context, _selected_parent_pool):
    _num_of_selected_parents = len(_selected_parent_pool)
    _mating_pool = []
    # print("-------------- START OF CROSSOVER AND MUTATION OF PARENT SELECTION --------------------------")
    for _index in range(_num_of_selected_parents):
        _index_1 = _index % _num_of_selected_parents
        _index_2 = (_index + 1) % _num_of_selected_parents
        _mating_pool.append(_selected_parent_pool[_index_1])
        _mating_pool.append(_selected_parent_pool[_index_2])

    for _chromosome in _mating_pool:
        _chromosome.print_chromosome("[GEN - " + str(_context.get_current_generation()) + "] :: Matting Parent :: ")

    _children_population_pool = []
    for _index in range(int(len(_mating_pool) / 2)):
        _parent_1 = _mating_pool[(2 * _index)]
        _parent_2 = _mating_pool[(2 * _index) + 1]
        _child_1, _child_2 = _parent_1.crossover(_parent_2)
        _child_1.mutation(_context.get_mutation_rate())
        _child_2.mutation(_context.get_mutation_rate())

        calculate_fitness_score(_context, _child_1)
        calculate_fitness_score(_context, _child_2)

        _children_population_pool.append(_child_1)
        _children_population_pool.append(_child_2)

    # print("-------------- AFTER CROSSOVER AND MUTATION AND BEFORE SORTING OF "
    #      "PARENT SELECTION --------------------------")

    # for _chromosome in _context.get_population_pool():
    #    _chromosome.print_chromosome("[GEN - " + str(_context.get_current_generation())
    #                                 + "] :: Before Sorting New Child :: ")

    # Sorting (in ascending order) all the chromosomes present in the population pool
    _context.initialize_population_pool(sorted(_children_population_pool, key=lambda element: element.fitness_score,
                                               reverse=True))
    print("-------------- AFTER CROSSOVER AND MUTATION AND SORTING OF PARENT SELECTION  --------------------------")

    for _chromosome in _context.get_population_pool():
        _chromosome.print_chromosome("[GEN - " + str(_context.get_current_generation()) + "] :: New Child :: ")

    return _context.get_population_pool()


def execute_generate_eight_queen_solution(_number_of_queens=8,
                                          _number_of_members_in_pool=10,
                                          _mutation_rate=0.1):
    _target_chromosome = Chromosome([0] * _number_of_queens)
    _is_valid = False
    _target_solution_found = False
    try:
        execution_context = GAContext(_number_of_members_in_pool, _number_of_queens, _mutation_rate)
        current_generation = 1
        _max_generation = 1000
        _initial_population_pool = generate_initial_population(execution_context)
        _max_scored_chromosome = _initial_population_pool[0]
        _max_possible_score = _number_of_queens * (_number_of_queens - 1)

        if _max_scored_chromosome.get_fitness_score() == _max_possible_score:
            _target_chromosome = _max_scored_chromosome
            _max_scored_chromosome.print_chromosome(
                "Generation #" + str(current_generation) + " Target Chromosome :: ")
            _target_solution_found = True
        else:

            while True:
                selected_parent_pool, max_score = run_parent_selection_process(execution_context)
                _initial_population_pool = execute_crossover_mutation(execution_context, selected_parent_pool)
                _max_scored_chromosome = _initial_population_pool[0]
                current_generation = current_generation + 1
                execution_context.set_current_generation(current_generation)

                if _max_scored_chromosome.get_fitness_score() == _max_possible_score:
                    _target_chromosome = _max_scored_chromosome
                    _max_scored_chromosome.print_chromosome(
                        "Generation #" + str(current_generation) + " Target Chromosome :: ")
                    _target_solution_found = True
                    break

                if current_generation > _max_generation:
                    print("Number of Maximum Generation limit crossed. Terminating the engine....")
                    _is_valid = False
                    break

                highest_score_candidate = _initial_population_pool[0]
                second_highest_score_candidate = _initial_population_pool[1]
                if current_generation % 10 == 0:
                    highest_score_candidate.print_chromosome("Generation #" + str(current_generation) + " ::: Highest ")
                    second_highest_score_candidate.print_chromosome(
                        "Generation #" + str(current_generation) + " ::: Second Highest ")

        if _target_solution_found:
            _target_chromosome.print_chromosome("Generation #" + str(current_generation) + " ::: Final ")
            _is_valid = validate_fitness_score(_target_chromosome.get_raw_data_sequence())
            if _is_valid:
                print("The final sequence {} is a valid solution.".format(_target_chromosome.get_raw_data_sequence()))
            else:
                print(
                    "The final sequence {} is a not valid solution.".format(_target_chromosome.get_raw_data_sequence()))

    except (TypeError, ValueError):
        raise

    return _target_chromosome, _target_solution_found


if __name__ == "__main__":

    try:
        _target_chromosome, is_target_solution_found = execute_generate_eight_queen_solution(8, 10, 0.1)

        if is_target_solution_found:
            print("The target sequence {}.".format(_target_chromosome.get_raw_data_sequence()))
        else:
            print("Unable to find any target sequence withing 2000 generations.")

    except (TypeError, ValueError):
        raise
