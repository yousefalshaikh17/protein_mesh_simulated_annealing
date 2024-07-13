import unittest
from tetmeshtools.scripts.temp.simulated_annealing import compute_probability_of_keeping

#@TODO work on naming
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
