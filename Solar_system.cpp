/// Please use the following command to compile the code:
// g++ -std=c++14 -o solar_system EE553 hw5_SolarSystem.cpp or just hit F5 for window VS code 
// use file solarsystem.dat to read the data
// Use regular C++ libraries and avoid using any non-standard libraries or functions.
// solorsystem.dat Mars Diam(m) changed from  6.792e6 to 6.792e7 due to error in the file , wrong output on person weight. 

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <random>
#include <sstream>
#include <thread> // Required for sleep functionality
#include <chrono> // Required for specifying time units
#include <conio.h> // For _getch()
using namespace std;

const double G = 6.674E-11;  // Gravitational constant
const double pi = 3.14159265358979323846;  // Pi constant

// Struct representing a 3D vector
struct Vec3d {
    double x, y, z;

    // Operator overloading for vector addition
    Vec3d operator+(const Vec3d& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }

    // Operator overloading for vector addition assignment
    Vec3d& operator+=(const Vec3d& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    // Operator overloading for scalar multiplication
    Vec3d operator*(double scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }

    // Overloading the stream insertion operator for Vec3d
    friend ostream& operator<<(ostream& os, const Vec3d& vec) {
        os << vec.x << ", " << vec.y << ", " << vec.z;
        return os;
    }
};

// Class representing a celestial body in the solar system
class Body {
private:
    string name;      // Name of the body
    string orbit;     // Orbit type (e.g., Sun)
    double mass;      // Mass of the body
    Vec3d pos;        // Position vector
    Vec3d v;          // Velocity vector
    Vec3d a;          // Acceleration vector
    double diameter;  // Diameter of the planet in meters
public:
    // Default constructor
    Body() : name("none"), orbit("none"), mass(0), pos({0, 0, 0}), v({0, 0, 0}), a({0, 0, 0}) {}

    // Parameterized constructor
    Body(string name, string orbit, double mass, Vec3d pos, Vec3d v, Vec3d a,double diam)
        : name(name), orbit(orbit), mass(mass), pos(pos), v(v), a(a),diameter(diam) {}

    // Getter and setter methods
    string getName() const { return name; }
    string getOrbit() const { return orbit; }
    double getMass() const { return mass; }
    Vec3d getPosition() const { return pos; }
    Vec3d getVelocity() const { return v; }
    Vec3d getAcceleration() const { return a; }
    double getDiameter() const { return diameter; }
    double getRadius() const { return diameter / 2; }  // Radius is half of the diameter

    void setName(const string& name) { this->name = name; }
    void setOrbit(const string& orbit) { this->orbit = orbit; }
    void setMass(double mass) { this->mass = mass; }
    void setPosition(const Vec3d& pos) { this->pos = pos; }
    void setVelocity(const Vec3d& v) { this->v = v; }
    void setAcceleration(const Vec3d& a) { this->a = a; }



    // Calculate gravity on the surface of the planet
    double calculateGravity() const {
        double radius = getRadius(); // Use radius instead of diameter
        if (radius <= 0) {
            cerr << "Error: Invalid radius for planet " << getName() << ".\n";
            return 0; // Avoid division by zero
        }
        if (mass <= 0) {
            cerr << "Error: Invalid mass for planet " << getName() << ".\n";
            return 0; // Avoid invalid mass
        }
        return G * mass / (radius * radius);
    }

    // Calculate a person's weight on this planet
    double calculateWeight(double personMass) const {
        double gravity = calculateGravity();
        if (gravity <= 0) {
            cerr << "Error: Could not calculate weight due to invalid gravity.\n";
            return 0; // Return 0 in case of error
        }
        return personMass * gravity; // Weight in Newtons
    }
  
      // Overloading the stream insertion operator for Body
    friend ostream& operator<<(ostream& os, const Body& body) {
        os << "Planet Name: " << body.name << "\n"
           << "  Orbit: " << body.orbit << "\n"
           << "  Mass: " << body.mass << " kg\n"
           << "  Position: (" << body.pos << ") m\n"
           << "  Velocity: (" << body.v << ") m/s\n"
           << "  Acceleration: (" << body.a << ") m/s²\n";
        return os;
    }
    // Static method to set acceleration for all bodies
    static void setAccelerations(vector<Body>& bodies, int acceleration) {
        Vec3d new_acceleration = {static_cast<double>(acceleration), static_cast<double>(acceleration), static_cast<double>(acceleration)};
        for (auto& body : bodies) {
            body.setAcceleration(new_acceleration);
        }
    }

    friend class SolarSystem;
};

// Class representing the solar system
class SolarSystem {
private:
    vector<Body> bodies;  // List of celestial bodies in the solar system
    double sunMass;       // Mass of the sun

public:
    // Constructor to initialize the solar system from a file
    SolarSystem(const string& filename) {
        ifstream solarfile(filename);
        if (!solarfile.is_open()) {
            cerr << "Error: File not found! Check if the path is correct.\n";
            return;
        }
    
        cout << "File is open and ready\n\n";
        string line;
        getline(solarfile, line); // Skip the header line
    
        // Variables to store body properties
        string name, orbit;
        double mass, diameter, x, y, z;
    
        // Seed the random number generator
        random_device rd;
        // 
        mt19937 gen(rd()); // random number generator based on the Mersenne Twister algorithm
        uniform_real_distribution<> posDis(-1e11, 1e11); // Position range: ±100 billion meters
        uniform_real_distribution<> angleDis(0, 2 * pi); // Angle range: 0 to 2π
    
        // Read each line from the file and create Body objects
        while (getline(solarfile, line)) {
            istringstream ss(line);
            ss >> name >> orbit >> mass >> diameter;

            
            // Assign random positions
            x = posDis(gen);
            y = posDis(gen);
            z = posDis(gen);
    
            if (name == "Sun") {
                sunMass = mass; // Store the Sun's mass
            }
    
            // Calculate velocity and acceleration for bodies orbiting the Sun
            double vx = 0, vy = 0, vz = 0, ax = 0, ay = 0, az = 0;
            if (orbit == "Sun") {
                double radius = sqrt(x * x + y * y + z * z); // Calculate distance from the Sun
                if (radius <= 0) {
                    cerr << "Error: Invalid radius for " << name << ".\n";
                    continue; // Skip this body
                }
    
                double orbitVelocity = sqrt(G * sunMass / radius); // Orbital velocity
                double centripetalAcceleration = (orbitVelocity * orbitVelocity) / radius; // Centripetal acceleration
    
                double angle = angleDis(gen); // Random orbital angle
    
                vx = orbitVelocity * cos(angle);
                vy = orbitVelocity * sin(angle);
                vz = 0;
    
                ax = centripetalAcceleration * cos(angle);
                ay = centripetalAcceleration * sin(angle);
                az = 0;
    
                // Print debug information
                cout << "================ ===================== =====================\n";
                cout << "Body name: " << name << " orbit: Sun\n"
                     << "Orbital velocity: " << orbitVelocity << " m/s\n"
                     << "Centripetal acceleration: " << centripetalAcceleration << " m/s²\n\n";
            }
    
            // Create a new Body object and add it to the list
            Body tempBody(name, orbit, mass, {x, y, z}, {vx, vy, vz}, {ax, ay, az}, diameter);
            bodies.push_back(tempBody);
        }
    
        // Close the file after reading
        solarfile.close();
    }
    
    // Method to update the solar system's state
    void stepForward(int acc) {
        Body::setAccelerations(bodies, acc);
    }

    // Overloading the stream insertion operator for SolarSystem
    friend ostream& operator<<(ostream& os, const SolarSystem& solarSystem) {
        for (const auto& body : solarSystem.bodies) {
            os << body << "\n";
        }
        return os;
    }
    // Method to update the solar system's state
    void update(double dt) {
        const double epsilon = 1e-3; // Minimum distance threshold to avoid singularities
    
        // Update accelerations due to gravitational forces
        for (auto& body : bodies) {
            Vec3d totalAcceleration = {0, 0, 0};
            for (const auto& other : bodies) {
                if (&body != &other) {
                    double dx = other.getPosition().x - body.getPosition().x;
                    double dy = other.getPosition().y - body.getPosition().y;
                    double dz = other.getPosition().z - body.getPosition().z;
    
                    double distance = sqrt(dx * dx + dy * dy + dz * dz);
                    if (distance < epsilon) {
                        cerr << "Warning: Very small distance between " << body.getName() 
                             << " and " << other.getName() << ". Skipping force calculation.\n";
                        continue;
                    }
    
                    // Gravitational force calculation
                    double force = (G * body.getMass() * other.getMass()) / (distance * distance);
                    totalAcceleration.x += force * (dx / distance) / body.getMass();
                    totalAcceleration.y += force * (dy / distance) / body.getMass();
                    totalAcceleration.z += force * (dz / distance) / body.getMass();
                }
            }
    
            // Update the body's acceleration
            body.setAcceleration(totalAcceleration);
        }
    
        // Update positions and velocities based on acceleration
        for (auto& body : bodies) {
            Vec3d newVelocity = body.getVelocity() + body.getAcceleration() * dt;
            Vec3d newPosition = body.getPosition() + newVelocity * dt;
    
            body.setVelocity(newVelocity);
            body.setPosition(newPosition);
        }
    }
    

     // Method to print the weights of a person on each planet
     void printWeights(double personMass) const {
        const double earthGravity = 9.81; // Earth's gravity in m/s²

        cout << "================ Your weight on Each Planet ================\n";
       
        for (const auto& body : bodies) {
            if (body.getOrbit() == "Sun") {  // Consider only planets orbiting the Sun
                double weight = body.calculateWeight(personMass);   // Calculate weight (N)
                double massEquivalent = weight / earthGravity;     // Calculate mass equivalent (kg)
                cout << "Planet: " << body.getName()
                 << ", Weight: " << weight << " N"
                 << ", Mass Equivalent: " << massEquivalent << " kg\n";
            }
        }
    }

    // Method to print the details of the solar system
    void printDetails() const {
        cout << "================ Solar System Details (First 3) ================\n\n";
        int count = 0;
        for (const auto& body : bodies) {
            if (count >= 3) break;  // Limit to first 3 bodies  to limit output for better visualization
           
            cout << "Planet Name    : " << body.getName() << "\n"
                 << "  Orbit        : " << body.getOrbit() << "\n"
                 << "  Mass         : " << body.getMass() << " kg\n"
                 << "  Diameter     : " << body.getDiameter() << " m\n"
                 << "  Position     : (" << body.getPosition() << ") m\n"
                 << "  Velocity     : (" << body.getVelocity() << ") m/s\n"
                 << "  Acceleration : (" << body.getAcceleration() << ") m/s²\n\n";
        
            count++; // Increment the counter    
        }
    }

    // Method to animate the printing of text
    void animatePrint(const std::string& text, int delayMs = 60) {
    for (char c : text) {
        std::cout << c << std::flush; // Print each character and flush the output
        std::this_thread::sleep_for(std::chrono::milliseconds(delayMs)); // Delay
    }
    std::cout << std::endl; // Move to the next line after printing
    }
    // Method to animate the printing of multiple lines
    void animatePrintLines(const std::vector<std::string>& lines, int delayMs = 500) {
        for (const std::string& line : lines) {
            std::cout << line << std::endl;
            std::this_thread::sleep_for(std::chrono::milliseconds(delayMs)); // Delay between lines
        }
    }
    // Method to get the list of bodies
    vector<Body>& getBodies() { return bodies; }
};


int main() {

    cout << "===========================================================\n";
    cout << "🌌 Welcome, Explorer of the Cosmos! 🌠\n";
    cout << "===========================================================\n\n";
    std::this_thread::sleep_for(std::chrono::seconds(2)); // Pause for 2 seconds

    cout << "Booting up the Solar System Analyzer... 🚀\n";
    std::this_thread::sleep_for(std::chrono::seconds(2)); // Pause for 2 seconds
    

    
    cout << "=====[ START ]=====\n";
     
    

    // Create a SolarSystem object from a file
    // SolarSystem s("solarsystem.dat");m
    // //
    
    
        // Initialize the SolarSystem object with the data file
    SolarSystem solarSystem("solarsystem.dat");  // Ensure the file path is correct

        // Print the details of the solar system
    string data2 =  "Printing Solar System Details:\n";
    solarSystem.animatePrint(data2);
    
    std::this_thread::sleep_for(std::chrono::seconds(2));


    // Print the details of the solar system
    solarSystem.printDetails();
    
    string data4 ="Changing Acceleration Value:\n";
    solarSystem.animatePrint(data4);
    int acc;
    cout << "Please enter the acceleration value in m/s²: ";
    cin >> acc;


    cout << "Hold on tight! Changing acceleration value...\n\n";
    std::this_thread::sleep_for(std::chrono::seconds(2));

    cout << "   Acceleration Value: " << acc << " m/s²\n\n";

    // Update the solar system's state with the new acceleration
    solarSystem.stepForward(acc);


    // Update the solar system's state
    double dt = 1.0;
    solarSystem.update(dt);

    cout << "Updated Solar System State:\n";

    string data3 =  "Printing New Solar System Details:\n";
    solarSystem.animatePrint(data3);
    std::this_thread::sleep_for(std::chrono::seconds(2));

    // Print the updated details of the solar system
    solarSystem.printDetails();



    cout << "===========================================================\n";

    //  cout << "🌌 Welcome to the Solar System Weight Calculator! 🌍\n";

    std::this_thread::sleep_for(std::chrono::seconds(2));

    string data = "🌌 Welcome to the Solar System Weight Calculator! 🌍\n";
    solarSystem.animatePrint(data);
    
    // Prompt user for their mass
    cout << "Please enter your mass in kilograms to see what you'd weigh on other planets: ";
    double personMass;
    cin >> personMass;
    
    
    cout << "\nHold on tight! Calculating your weight across the Solar System... 🚀\n\n";
    std::this_thread::sleep_for(std::chrono::seconds(2));
    solarSystem.printWeights(personMass);

    solarSystem.animatePrintLines({
        "===========================================================",
        "🌌 Thank you for using the Solar System Weight Calculator! 🌌",
        "We hope you enjoyed your journey through the cosmos! 🚀",
        "==========================================================="
    }, 1000);

    
    cout << "\n\n";


    
 


    cout << "=====================[ Mission Complete! ]=================\n";
  


    return 0;
}
