import unittest
from tetmeshtools.scripts.simulated_annealing.simulated_annealing import compute_probability_of_keeping, mutate
from tetmeshtools.meshtools.tetgenstructs import NodePoint

def is_within_threshold(original, mutated, threshold):
    diff_x = abs(original.x - mutated.x)
    diff_y = abs(original.y - mutated.y)
    diff_z = abs(original.z - mutated.z)
    return diff_x <= threshold and diff_y <= threshold and diff_z <= threshold

class TestSimulatedAnnealing(unittest.TestCase):
    """
    tests for simulated annealing
    """

    def setUp(self):
        """
        build a full test class
        """
        pass

    def tearDown(self):
        """
        clean up
        """
        pass

    def test_probability_of_keeping1(self):
        """
        test probability of keeping at a lower temperature
        """
        old_fitness = 30
        new_fitness = 40
        temperature = 80
        k = 0.15
        expected_result = 0.435

        output = compute_probability_of_keeping(old_fitness, new_fitness, temperature, k=k)
        self.assertAlmostEqual(output, expected_result, delta=0.001, msg="simulated_annealing.compute_probability_of_keeping")

    def test_probability_of_keeping2(self):
        """
        test probability of keeping at a lower temperature
        """
        old_fitness = 30
        new_fitness = 40
        temperature = 10
        k = 0.15
        expected_result = 0.001 # Should be drastically lower due to lower temperature.

        output = compute_probability_of_keeping(old_fitness, new_fitness, temperature, k=k)
        self.assertAlmostEqual(output, expected_result, delta=0.001, msg="simulated_annealing.compute_probability_of_keeping")

    def test_probability_of_keeping3(self):
        """
        test probability of keeping at a lower temperature with lower delta fitness
        """
        old_fitness = 23
        new_fitness = 25
        temperature = 10
        k = 0.15
        expected_result = 0.263 # Smaller delta fitness should have higher probability

        output = compute_probability_of_keeping(old_fitness, new_fitness, temperature, k=k)
        self.assertAlmostEqual(output, expected_result, delta=0.001, msg="simulated_annealing.compute_probability_of_keeping")


    def test_mutate1(self):
        """
        test mutation of node at a low temperature
        """
        temperature = 5
        mutation_multiplier = 1
        original_node = NodePoint(1, 34.2, 43.6, 13.1)
        for _ in range(10):
            mutated_node = mutate(original_node, temperature, mutation_multiplier=mutation_multiplier)
            # If mutated node is identical, try once more.
            if mutated_node == original_node:
                mutated_node = mutate(original_node, temperature, mutation_multiplier=mutation_multiplier)
                self.assertNotEqual(original_node, mutated_node, msg="simulated_annealing.test_mutate (Node was not mutated.)")
            self.assertTrue(is_within_threshold(original_node, mutated_node, temperature * mutation_multiplier), msg="simulated_annealing.test_mutate (Mutation went over the limits)")

    def test_mutate2(self):
        """
        test mutation of node at a higher temperature
        """
        temperature = 30
        mutation_multiplier = 0.7
        original_node = NodePoint(1, 662.2, -139.6, -439.1)
        for _ in range(10):
            mutated_node = mutate(original_node, temperature, mutation_multiplier=mutation_multiplier)
            # If mutated node is identical, try once more.
            if mutated_node == original_node:
                mutated_node = mutate(original_node, temperature, mutation_multiplier=mutation_multiplier)
                self.assertNotEqual(original_node, mutated_node, msg="simulated_annealing.test_mutate (Node was not mutated.)")
            self.assertTrue(is_within_threshold(original_node, mutated_node, temperature * mutation_multiplier), msg="simulated_annealing.test_mutate (Mutation went over the limits)")
