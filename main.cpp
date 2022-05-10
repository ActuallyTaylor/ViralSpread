#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <random>
#include <cmath>
#include<algorithm>
#include <execution>

/*
 An Implementation of a variation of the SIERS model for infectious disease.
 This simulation attempts to model the spread, infection, and evolution of viral diseases across a population (The World)
 */

using namespace std;
int rollingPersonID { 0 };
int rollingVirusID { 0 };

/*
 * World Constants
 */
// The number of daily social interactions. Based on resarch from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6113687/
int numberOfInteractionsPerTrip = 3;
// The average number of times somebody travels. Data calculated on findings. https://www.bts.gov/archive/publications/highlights_of_the_2001_national_household_travel_survey/section_02
// & https://www.researchgate.net/publication/279853330_Extended_Range_Electric_Vehicle_Driving_and_Charging_Behavior_Observed_Early_in_the_EV_Project
// Trips are defined as moving from one address to another. This means, they usually come in multiples of 2.
int numberOfDailyTrips = 4;
int averageTripDistanceMiles = 8;

float randomFloat() {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<float> dist(0.0, 1.0);
    return dist(mt);
}

int randomInt(int min, int max) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<> distr(min, max);
    return distr(mt);
}

int randomWeightedInt(initializer_list<double> weights,  initializer_list<int> numbers) {
    random_device rd;
    mt19937 gen(rd());
    discrete_distribution<> d(weights);
    map<int, int> map;
    int index = d(gen);
    return numbers.begin()[index];
}

template <typename T>
int indexOf(T object, vector<T> vector) {
    for(int n = 0; n < vector.capacity(); n ++) {
        if(vector[n] == object) {
            return n;
        }
    }
    return -1;
}

// https://stackoverflow.com/a/9345144/14886210
template<class BidiIter>
BidiIter random_unique(BidiIter begin, BidiIter end, size_t num_random) {
    size_t left = distance(begin, end);
    while (num_random--) {
        BidiIter r = begin;
        advance(r, rand()%left);
        swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}

enum VirusType {
    Lysogenic, // Long-lasting, can stay dormant for a while
    Lytic // Short, usually takes effect right away
};

enum InfectionMethod {
    Touch,
    Saliva,
    Coughing,
    Sneezing,
    SexualContact,
    Contamination,
    Insects
};

enum DiseaseStage {
    Susceptible, // Can become infected
    Infected, // Infected, but unable to infect others
    Infectious, // Able to infect others
    Recovered, // Immune to the disease
    Immune // A temporary immunity, comes from birth, immunity decreases after birth
};

struct Virus {
    string name;
    int id;
    vector<InfectionMethod> infectionMethods { };

    float infectionRate { 0.025 }; // A percentage out of 100
    int latencyPeriod { 7 }; // The period at which it takes to cross from infected -> infectious

    Virus(string name, vector<InfectionMethod> infectionMethods) {
        this->name = std::move(name);
        this->id = rollingVirusID;
        this->infectionRate = randomFloat();
        this->infectionMethods = std::move(infectionMethods);
        rollingVirusID ++;
    }
};

struct City;
struct CityMile;

struct Person {
    int id;
    int x { 0 };
    int y { 0 };
    int parentChunkIndex;

    int tripsToday { 0 };
    int tripsCounter { 0 };
    optional<Virus> infectiousVirus;

    map<int, float> strainResistance { }; // Virus ID | Resistance percent 0-1

    // https://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    // This function will randomly select a predisposed health condition based on a weighted random number generator.
    static int randomHealthCondition() {
        random_device rd;
        mt19937 gen(rd());
        discrete_distribution<> d({10, 15, 30, 30, 15});
        map<int, int> map;
        return d(gen);
    }

    DiseaseStage stageOfInfection { Susceptible };
    int StartingHealthCondition = randomHealthCondition();

    float currentHealthCondition { 1 }; // Start off with perfect health, 0 - 1.0
    float chanceOfInfection() { // Need to actually implement, 0 - 1.0
        return currentHealthCondition / 1;
    }

    explicit Person(int id, int parentChunkIndex) {
        this->id = id;
        this->parentChunkIndex = parentChunkIndex;
    }

    inline bool operator == (const Person& rhs) const {
        return id == rhs.id;
    }

    bool infect(Virus virus, bool override = false) {
        bool canInfect = chanceOfInfection() <= virus.infectionRate;
        if (canInfect || override) {
            infectiousVirus = virus;
        }
        return canInfect;
    }
};

struct CityMile {
    int id;
    int xPosition;
    int yPosition;
    int neighborIndices [8];
    vector<Person*> people { };

    CityMile() {
        id = 0;
        xPosition = 0;
        yPosition = 0;
    }

    CityMile(int id, int xPosition, int yPosition) {
        this->id = id;
        this->xPosition = xPosition;
        this->yPosition = yPosition;
    }

    int randomNeighborIndex() {
        vector<int> possibleIndices;

        for(int x = 0; x < 8; x++) {
            if (neighborIndices[x] != -1) {
                possibleIndices.push_back(x);
            }
        }

        if(!possibleIndices.empty()) {
            int random = randomInt(0, possibleIndices.size() - 1);
            return possibleIndices[random];
        } else {
            return 0;
        }
    }
};

struct City {
    // Function to check if a number is a perfect square or not
    static bool isPerfect(int number) {
        return (sqrt(number) - floor(sqrt(number))) == 0;
    }

    // Adapted from https://www.geeksforgeeks.org/closest-perfect-square-and-its-distance/
    // Function to find the closest perfect square
    static int getClosestPerfectSquare(int number) {
        // Check if it is already a perfect number
        if (isPerfect(number)) {
            return number;
        }

        // Variables to store the perfect squares
        int squareAbove = -1, squareBelow = -1;
        int counter;

        // Find first perfect square above the given number
        counter = number + 1;
        while (true) {
            if (isPerfect(counter)) {
                squareAbove = counter;
                break;
            } else {
                counter++;
            }
        }

        // Find first perfect square below the given number
        counter = number - 1;
        while (true) {
            if (isPerfect(counter)) {
                squareBelow = counter;
                break;
            } else {
                counter--;
            }
        }

        // Return the closest square
        if ((squareAbove - number) > (number - squareBelow)) {
            return squareBelow;
        } else {
            return squareAbove;
        }
    }

    string name;
    int population { };
    float populationDensity { };
    float squareMileage { };

    int arraySideLength { 0 };
    vector<CityMile> cityChunks { };
    vector<Person> people { };

    City() {
        this->name = "";
        this->population = 0;
    }

    City(string name, int population, float populationDensity, float squareMileage) {
        this->name = std::move(name);
        this->population = population;
        this->populationDensity = populationDensity;
        this->squareMileage = squareMileage;

        int roundedPopulationDensity = ceil(populationDensity);
        int numberOfChunks = ceil(squareMileage);

        arraySideLength = sqrt(getClosestPerfectSquare(numberOfChunks));
        int x { 0 };
        int y { 0 };

        // Create a chunk for every needed chunk.
        for(int index = 0; index < numberOfChunks; index++) {
            for(int pIndex = 0; pIndex < roundedPopulationDensity; pIndex++) {
                Person tempPerson = Person(rollingPersonID, index);
                people.push_back(tempPerson);
                rollingPersonID ++;
            }
            cityChunks.emplace_back(index, x, y);

            if (x >= arraySideLength - 1) {
                x = 0;
                y ++;
            } else {
                x ++;
            }

            // Make each chunk aware of its neighbors
            // Check to see if we are on the edge of the square
            if(index % arraySideLength == 1) {
                // We are on the left side
                // North Western Chunk
                cityChunks[index].neighborIndices[0] = -1;

                // Northern Chunk
                if (index - arraySideLength >= 0) {
                    cityChunks[index].neighborIndices[1] = index - arraySideLength;
                } else {
                    cityChunks[index].neighborIndices[1] = -1;
                }

                // North Eastern Chunk
                if (index - arraySideLength + 1 >= 0) {
                    cityChunks[index].neighborIndices[2] = index - arraySideLength + 1;
                } else {
                    cityChunks[index].neighborIndices[2] = -1;
                }

                // Western Chunk
                cityChunks[index].neighborIndices[3] = -1;

                // Eastern Chunk
                if (index + 1 >= 0 && index + 1 <= cityChunks.size()) {
                    cityChunks[index].neighborIndices[4] = index + 1;
                } else {
                    cityChunks[index].neighborIndices[4] = -1;
                }

                // South Western Chunk
                cityChunks[index].neighborIndices[5] = -1;

                // Southern Chunk
                if(index + arraySideLength <= cityChunks.size()) {
                    cityChunks[index].neighborIndices[6] = index + arraySideLength;
                } else {
                    cityChunks[index].neighborIndices[6] = -1;
                }

                // South Eastern Chunk
                if(index + arraySideLength + 1 <= cityChunks.size()) {
                    cityChunks[index].neighborIndices[7] = index + arraySideLength + 1;
                } else {
                    cityChunks[index].neighborIndices[7] = -1;
                }
            } else if (index % arraySideLength == 0) {
                // We are on the right side
                // North Western Chunk
                if (index - arraySideLength - 1 >= 0) {
                    cityChunks[index].neighborIndices[0] = index - arraySideLength - 1;
                } else {
                    cityChunks[index].neighborIndices[0] = -1;
                }

                // Northern Chunk
                if (index - arraySideLength >= 0) {
                    cityChunks[index].neighborIndices[1] = index - arraySideLength;
                } else {
                    cityChunks[index].neighborIndices[1] = -1;
                }

                // North Eastern Chunk
                cityChunks[index].neighborIndices[2] = -1;

                // Western Chunk
                if (index - 1 >= 0) {
                    cityChunks[index].neighborIndices[3] = index - 1;
                } else {
                    cityChunks[index].neighborIndices[3] = -1;
                }

                // Eastern Chunk
                cityChunks[index].neighborIndices[4] = -1;

                // South Western Chunk
                if (index + arraySideLength - 1 >= 0) {
                    cityChunks[index].neighborIndices[5] = index + arraySideLength - 1;
                } else {
                    cityChunks[index].neighborIndices[5] = -1;
                }

                // Southern Chunk
                if(index + arraySideLength <= cityChunks.size()) {
                    cityChunks[index].neighborIndices[6] = index + arraySideLength;
                } else {
                    cityChunks[index].neighborIndices[6] = -1;
                }

                // South Eastern Chunk
                cityChunks[index].neighborIndices[7] = -1;
            } else {
                // We are in the middle somewhere
                // North Western Chunk
                if (index - arraySideLength - 1 >= 0) {
                    cityChunks[index].neighborIndices[0] = index - arraySideLength - 1;
                } else {
                    cityChunks[index].neighborIndices[0] = -1;
                }

                // Northern Chunk
                if (index - arraySideLength >= 0) {
                    cityChunks[index].neighborIndices[1] = index - arraySideLength;
                } else {
                    cityChunks[index].neighborIndices[1] = -1;
                }

                // North Eastern Chunk
                if (index - arraySideLength + 1 >= 0) {
                    cityChunks[index].neighborIndices[2] = index - arraySideLength + 1;
                } else {
                    cityChunks[index].neighborIndices[2] = -1;
                }

                // Western Chunk
                if (index - 1 >= 0) {
                    cityChunks[index].neighborIndices[3] = index - 1;
                } else {
                    cityChunks[index].neighborIndices[3] = -1;
                }

                // Eastern Chunk
                if (index + 1 >= 0 && index + 1 <= cityChunks.size()) {
                    cityChunks[index].neighborIndices[4] = index + 1;
                } else {
                    cityChunks[index].neighborIndices[4] = -1;
                }

                // South Western Chunk
                if (index + arraySideLength - 1 >= 0) {
                    cityChunks[index].neighborIndices[5] = index + arraySideLength - 1;
                } else {
                    cityChunks[index].neighborIndices[5] = -1;
                }

                // Southern Chunk
                if(index + arraySideLength <= cityChunks.size()) {
                    cityChunks[index].neighborIndices[6] = index + arraySideLength;
                } else {
                    cityChunks[index].neighborIndices[6] = -1;
                }

                // South Eastern Chunk
                if(index + arraySideLength + 1 <= cityChunks.size()) {
                    cityChunks[index].neighborIndices[7] = index + arraySideLength + 1;
                } else {
                    cityChunks[index].neighborIndices[7] = -1;
                }
            }
        }
    }
};

struct World {
//    vector<City> cities { };
    City cities [1];
    vector<Virus> viruses { };

    void infect(City& city, Virus& virus) {
        city.people[randomInt(0, city.people.size())].infect(virus, true);
    }

    void simulate() {
        for(City& city: cities) {
            // Loop over all the possible trip numbers
            int tripCount = 9;
            for(int trip = 0; trip < tripCount; trip++) {
                // Simulate disease spread in parallel threads.
                #pragma omp parallel for
                for(auto chunkIterator = begin (city.cityChunks); chunkIterator != end (city.cityChunks); ++chunkIterator) {
                    for(auto personIterator = begin (chunkIterator->people); personIterator != end (chunkIterator->people); ++personIterator) {
                        int numberOfInteractions  =  numberOfInteractionsPerTrip * (*personIterator)->tripsToday;
                        if (numberOfInteractions > chunkIterator->people.size()) {
                            numberOfInteractions = chunkIterator->people.size();
                        }

                        vector<Person*> interactions = chunkIterator->people;
                        random_unique(interactions.begin(), interactions.end(), numberOfInteractions);

                        for(int n = 0; n < numberOfInteractions; n++) {
                            Person* possibleInfector = interactions[n];
                            // If we are looking at ourselves, skip.
                            if (possibleInfector->id == (*personIterator)->id) {
                                continue;
                            }
                            // Check to see if the other person is infectious.
                            if(possibleInfector->stageOfInfection == DiseaseStage::Infectious) {
                                if(possibleInfector->infectiousVirus.has_value()) {
                                    (*personIterator)->infect(possibleInfector->infectiousVirus.value());
                                    cout << "Attempt to infect" << endl;
                                }
                            }
                        }
                    }

                    // Clear people vector to allow for movement to happen
                    chunkIterator->people.clear();
                }

                // Move person to a new chunk
                if (trip != tripCount - 1) {
                    for(auto personIterator = begin (city.people); personIterator != end (city.people); ++personIterator) {
                        if(trip == 0) {
                            int dailyTrips = numberOfDailyTrips + randomWeightedInt({5, 15, 60, 15, 5},
                                                                                    {-4, -2, 0, 2, 4});
                            personIterator->tripsToday = dailyTrips;
                            personIterator->tripsCounter = dailyTrips;
                        }

                        if (personIterator->tripsCounter < trip) {
                            int travelDistance = averageTripDistanceMiles +
                                                 randomWeightedInt({2.5, 10, 20, 30, 20, 10, 5, 2.5},
                                                                   {-5, -3, -1, 0, 1, 3, 5, 10});
                            for (int x = 0; x < travelDistance; x++) {
                                int parentChunkIndex = personIterator->parentChunkIndex;
                                int randomNeighborIndex = city.cityChunks[parentChunkIndex].randomNeighborIndex();
                                int indexInParentChunk = indexOf((*personIterator), city.people);
                                personIterator->parentChunkIndex = randomNeighborIndex;
                                personIterator->tripsCounter --;
                            }
                        }

                        city.cityChunks[personIterator->parentChunkIndex].people.push_back(&(*personIterator));
                    }
                }
            }
        }
    }
};

int main() {
    srand (static_cast <unsigned> (time(nullptr)));
    World Earth = World();

    // Read in cities from data (Currently Just One City)
    Earth.cities[0] = City("Philadelphia", 1000, 10, 10);
    Earth.viruses.emplace_back(Virus("Argo 1", {Touch, Saliva, Coughing, Sneezing, SexualContact, Contamination, Insects}));
    Earth.infect(Earth.cities[0], Earth.viruses[0]);

    Earth.simulate();

    return 0;
};

/* Check to see if weighted numbers are working
    double numb0 { 0 };
    double numb1 { 0 };
    double numb2 { 0 };
    double numb3 { 0 };
    double numb4 { 0 };

    int count = 1000;
    for(int x = 0; x < count; x++) {
        int condition = randomWeightedInt({15, 30, 30, 15, 10}, {0, 1, 2, 3, 4});
        if (condition == 0) {
            numb0 ++;
        } else if (condition == 1) {
            numb1 ++;
        } else if (condition == 2) {
            numb2 ++;
        } else if (condition == 3) {
            numb3 ++;
        } else if (condition == 4) {
            numb4 ++;
        }
    }

    cout << "0: " << numb0 << " " << ((numb0 / double(count)) * 100) << "%" << endl;
    cout << "1: " << numb1 << " " << ((numb1 / double(count)) * 100) << "%" << endl;
    cout << "2: " << numb2 << " " << ((numb2 / double(count)) * 100) << "%" << endl;
    cout << "3: " << numb3 << " " << ((numb3 / double(count)) * 100) << "%" << endl;
    cout << "4: " << numb4 << " " << ((numb4 / double(count)) * 100) << "%" << endl;

    OUTPUT:
    0: 145 14.5%
    1: 316 31.6%
    2: 302 30.2%
    3: 128 12.8%
    4: 109 10.9%
 */
