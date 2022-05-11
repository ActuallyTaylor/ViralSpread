#include <iostream>
#include <fstream>
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
double roundTripChance = 0.6;

// https://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
// These functions will randomly select a random numbers.
double randomDouble(double min, double max) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(min, max);
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
    int index = d(gen);
    return numbers.begin()[index];
}

double randomWeightedDouble(initializer_list<double> weights,  initializer_list<double> numbers) {
    random_device rd;
    mt19937 gen(rd());
    discrete_distribution<> d(weights);
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

    double infectionRate { 0.025 }; // A percentage out of 100
    int latencyPeriod { 3 }; // The period at which it takes to cross from infected -> infectious
    int infectionPeriod { 14 }; // The length of time that a person is infected for
    int recoveryPeriod { 4 }; // The length of the recovery period

    Virus(string name, vector<InfectionMethod> infectionMethods) {
        this->name = std::move(name);
        this->id = rollingVirusID;
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

    // Dealing with moving chunks
    int tripStartChunkIndex { 0 };
    int tripLength { 0 };
    int roundTripCounter { 0 };
    bool isOnRoundTrip { false };

    int tripsToday { 0 };
    int tripsCounter { 0 };
    optional<Virus> infectiousVirus;

    map<int, double> strainResistance { }; // Virus ID | Resistance percent 0-1

    DiseaseStage stageOfInfection { Susceptible };
    int stageCounter { 0 };

    double currentHealthCondition { 1.0 };

    explicit Person(int id, int parentChunkIndex) {
        this->id = id;
        this->parentChunkIndex = parentChunkIndex;
        currentHealthCondition = randomWeightedDouble({0.1, 0.15, 0.3, 0.3, 0.1}, {0.1, 0.4, 0.5, 0.8, 1.0});
    }

    inline bool operator == (const Person& rhs) const {
        return id == rhs.id;
    }

    // Advances the infection through the stages fo the SIERS model
    void advanceInfection() {
        if (infectiousVirus.has_value()) {
            if(stageOfInfection == Infected) {
                if(infectiousVirus.value().latencyPeriod < stageCounter) {
                    stageOfInfection = Infectious;
                    stageCounter = 0;
                }
            } else if (stageOfInfection == Infectious) {
                if(infectiousVirus.value().infectionPeriod < stageCounter) {
                    stageOfInfection = Recovered;
                    stageCounter = 0;
                }
            } else if(stageOfInfection == Recovered) {
                if(infectiousVirus.value().recoveryPeriod < stageCounter) {
                    stageOfInfection = Susceptible;
                    stageCounter = 0;
                    infectiousVirus.reset();
                }
            }
            stageCounter ++;
        }
    }

    // Calcualtes the chance of infection. Based on health condition and rate of transfer of the virus. Along with the type of interaction.
    double chanceOfInfection() { // Need to actually implement, 0 - 1.0
        double chanceToInfect = currentHealthCondition / 20;
        return chanceToInfect;
    }

    bool infect(Virus virus, bool override = false) {
        bool canInfect = chanceOfInfection() <= virus.infectionRate;
        if(stageOfInfection == Recovered || stageOfInfection == Infectious || stageOfInfection == Infected) {
            canInfect = false;
        }
        if (canInfect || override) {
            infectiousVirus = virus;
            stageOfInfection = Infected;
        }
        return canInfect;
    }

    bool isInfected() {
        return infectiousVirus.has_value();
    }
};

struct CityMile {
    int id;
    int xPosition;
    int yPosition;
    int neighborIndices [8] = { -1, -1, -1, -1, -1, -1, -1, -1 };
    vector<Person*> people { };
    vector<int> possibleIndices;

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

    void setPossibleIndicies() {
        for(int x = 0; x < 8; x++) {
            if (neighborIndices[x] != -1) {
                possibleIndices.push_back(x);
            }
        }
    }

    // Picks a random index in the list of valid neighbors.
    int randomNeighborIndex() {
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
    double populationDensity { };
    double squareMileage { };

    int arraySideLength { 0 };
    vector<CityMile> cityChunks { };
    vector<Person> people { };

    City() {
        this->name = "";
        this->population = 0;
    }

    City(string name, int population, double populationDensity, double squareMileage) {
        this->name = std::move(name);
        this->population = population;
        this->populationDensity = populationDensity;
        this->squareMileage = squareMileage;

        int roundedPopulationDensity = ceil(populationDensity);
        int numberOfChunks = ceil(squareMileage);
        int arraySize = numberOfChunks - 1;

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

            // Make each chunk aware of its neighbors
            // Check to see if we are on the edge of the square
            if(x == 0) { // Left Edge
                if(index - arraySideLength < 0) {
//                    cout << "Chunk Index: " << index << " Top Left Edge" << endl;
                    setIndex(index, 4, index + 1, arraySize);
                    setIndex(index, 6, index + arraySideLength, arraySize);
                    setIndex(index, 7, (index + arraySideLength) + 1, arraySize);
                } else if (index + arraySideLength > numberOfChunks - 1) {
//                    cout << "Chunk Index: " << index << " Bottom Left Edge" << endl;
                    setIndex(index, 1, index - arraySideLength, arraySize);
                    setIndex(index, 2, (index - arraySideLength) + 1, arraySize);
                    setIndex(index, 4, index + 1, arraySize);
                } else {
//                    cout << "Chunk Index: " << index << " Left Edge" << endl;
                    setIndex(index, 1, index - arraySideLength, arraySize);
                    setIndex(index, 2, (index - arraySideLength) + 1, arraySize);
                    setIndex(index, 4, index + 1, arraySize);
                    setIndex(index, 6, index + arraySideLength, arraySize);
                    setIndex(index, 7, (index + arraySideLength) + 1, arraySize);
                }
            } else if (x == arraySideLength - 1) { // Right Edge
                if(index - arraySideLength < 0) {
//                    cout << "Chunk Index: " << index << " Top Right Edge" << endl;
                    setIndex(index, 3, index - 1, arraySize);
                    setIndex(index, 5, (index + arraySideLength) - 1, arraySize);
                    setIndex(index, 6, index + arraySideLength, arraySize);
                } else if (index + arraySideLength > arraySize) {
//                    cout << "Chunk Index: " << index << " Bottom Right Edge" << endl;
                    setIndex(index, 0, (index - arraySideLength) - 1, arraySize);
                    setIndex(index, 1, index - arraySideLength, arraySize);
                    setIndex(index, 3, index - 1, arraySize);
                } else {
//                    cout << "Chunk Index: " << index << " Right Edge" << endl;
                    setIndex(index, 0, (index - arraySideLength) - 1, arraySize);
                    setIndex(index, 1, index - arraySideLength, arraySize);
                    setIndex(index, 3, index - 1, arraySize);
                    setIndex(index, 5, (index + arraySideLength) - 1, arraySize);
                    setIndex(index, 6, index + arraySideLength, arraySize);
                }
            } else { // No Edge
                if(index - arraySideLength < 0) {
//                    cout << "Chunk Index: " << index << " Top Edge" << endl;
                    setIndex(index, 3, index - 1, arraySize);
                    setIndex(index, 4, index + 1, arraySize);
                    setIndex(index, 5, (index + arraySideLength) - 1, arraySize);
                    setIndex(index, 6, index + arraySideLength, arraySize);
                    setIndex(index, 7, (index + arraySideLength) + 1, arraySize);
                } else if (index + arraySideLength > arraySize) {
//                    cout << "Chunk Index: " << index << " Bottom Edge" << endl;
                    setIndex(index, 0, (index - arraySideLength) - 1, arraySize);
                    setIndex(index, 1, index - arraySideLength, arraySize);
                    setIndex(index, 2, (index - arraySideLength) + 1, arraySize);
                    setIndex(index, 3, index - 1, arraySize);
                    setIndex(index, 4, index + 1, arraySize);
                } else {
//                    cout << "Chunk Index: " << index << " No Edge" << endl;
                    setIndex(index, 0, (index - arraySideLength) - 1, arraySize);
                    setIndex(index, 1, index - arraySideLength, arraySize);
                    setIndex(index, 2, (index - arraySideLength) + 1, arraySize);
                    setIndex(index, 3, index - 1, arraySize);
                    setIndex(index, 4, index + 1, arraySize);
                    setIndex(index, 5, (index + arraySideLength) - 1, arraySize);
                    setIndex(index, 6, index + arraySideLength, arraySize);
                    setIndex(index, 7, (index + arraySideLength) + 1, arraySize);
                }
            }

            // Register all the possible indices for the city
            cityChunks[index].setPossibleIndicies();

            if (x >= arraySideLength - 1) {
//                cout << index << endl;
                x = 0;
                y ++;
            } else {
//                cout << index;
                x ++;
            }
        }
    }

    void setIndex(int cityIndex, int arrayIndex, int chunkIndex, int arraySize) {
        if(chunkIndex < arraySize && chunkIndex >= 0) {
            cityChunks[cityIndex].neighborIndices[arrayIndex] = chunkIndex;
        }
    }
};

struct World {
//    vector<City> cities { };
    City cities [1];
    vector<Virus> viruses { };

    void infect(City& city, Virus& virus) {
        int randIndex = randomInt(0, city.people.size());
        city.people[randIndex].infect(virus, true);
    }

    void simulate(int dayNumb) {
        for(City& city: cities) {
            // Loop over all the possible trip numbers
            int tripCount = 9;
            for (auto personIterator = begin (city.people); personIterator != end (city.people); ++personIterator) {
                (*personIterator).advanceInfection();
            }

            string csvString = "";
            for(int trip = 0; trip < tripCount; trip++) {
                // Simulate disease spread in parallel threads.
                #pragma omp parallel for
                for(auto chunkIterator = begin (city.cityChunks); chunkIterator != end (city.cityChunks); ++chunkIterator) {
                    int numbSusceptible = 0;
                    int numbInfected = 0;
                    int numbInfectious = 0;
                    int numbRecovered = 0;
                    int numbImmune = 0;

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
                                }
                            }
                        }

                        DiseaseStage stage = (*personIterator)->stageOfInfection;
                        if(stage == Susceptible) {
                            numbSusceptible ++;
                        } else if (stage == Infected) {
                            numbInfected ++;
                        } else if (stage == Infectious) {
                            numbInfectious ++;
                        } else if (stage == Recovered) {
                            numbRecovered ++;
                        } else if (stage == Immune) {
                            numbImmune ++;
                        }
                    }

                    // Clear people vector to allow for movement to happen, we only want to clear if we are about to move people
                    if (trip != tripCount - 1) {
                        chunkIterator->people.clear();
                    } else {
                        // + "," + to_string(trip)
                        csvString += to_string(dayNumb) + "," +city.name + "," + to_string(city.people.size()) + "," + to_string(chunkIterator->id) + "," + to_string(chunkIterator->people.size()) + "," + to_string(numbSusceptible) + "," + to_string(numbInfected) + "," + to_string(numbInfectious) + "," + to_string(numbRecovered) + "," + to_string(numbImmune) + "\n";
                    }
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
                            if (!personIterator->isOnRoundTrip) { // We are not on a return trip so perform a move normally
                                if(randomDouble(0, 1) <= roundTripChance) {
                                    personIterator->isOnRoundTrip = true;
                                    personIterator->roundTripCounter = 0;
                                    personIterator->tripLength = randomWeightedInt({0.5, 0.3, 0.1},{2, 4, 6}); // Pick the number length of the trip based on trip cycles.
                                    personIterator->tripStartChunkIndex = personIterator->parentChunkIndex;
                                }

                                int travelDistance = averageTripDistanceMiles +
                                                     randomWeightedInt({2.5, 10, 20, 30, 20, 10, 5, 2.5},
                                                                       {-5, -3, -1, 0, 1, 3, 5, 10});
                                for (int x = 0; x < travelDistance; x++) {
                                    int parentChunkIndex = personIterator->parentChunkIndex;
                                    int randomNeighborIndex = city.cityChunks[parentChunkIndex].randomNeighborIndex();
                                    personIterator->parentChunkIndex = randomNeighborIndex;
                                    personIterator->tripsCounter --;
                                }
                                personIterator->roundTripCounter ++;
                            } else { // We are on a round trip
                                if(personIterator->roundTripCounter >= personIterator->tripLength - 1) { // We need to be home on our next trip
                                    personIterator->parentChunkIndex = personIterator->tripStartChunkIndex;
                                    personIterator->tripStartChunkIndex = 0;
                                    personIterator->tripLength = 0;
                                    personIterator->isOnRoundTrip = false;

                                    personIterator->tripsCounter --;
                                } else { // We can keep moving
                                    int travelDistance = averageTripDistanceMiles +
                                                         randomWeightedInt({2.5, 10, 20, 30, 20, 10, 5, 2.5},
                                                                           {-5, -3, -1, 0, 1, 3, 5, 10});
                                    for (int x = 0; x < travelDistance; x++) {
                                        int parentChunkIndex = personIterator->parentChunkIndex;
                                        int randomNeighborIndex = city.cityChunks[parentChunkIndex].randomNeighborIndex();
                                        personIterator->parentChunkIndex = randomNeighborIndex;
                                        personIterator->tripsCounter --;
                                    }
                                    personIterator->roundTripCounter ++;
                                }
                            }
                        }

                        city.cityChunks[personIterator->parentChunkIndex].people.push_back(&(*personIterator));
                    }
                }
            }
            ofstream dataCSV;
            dataCSV.open ("./Daily_Infection_Numbers.csv", ios_base::app);
            dataCSV << csvString;
            dataCSV.close();
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

    ofstream dataCSV;
    dataCSV.open ("./Daily_Infection_Numbers.csv");
    //Trip,
    dataCSV << "Date,City,City Population,Chunk Number,Chunk Population,Susceptible,Infected,Infectious,Recovered,Immune\n";
    dataCSV.close();

    int daysToSimulate { 365 };
    for(int x = 0; x < daysToSimulate; x++) {
        Earth.simulate(x);
        cout << "Finished day: " << x << endl;
    }

    return 0;
};