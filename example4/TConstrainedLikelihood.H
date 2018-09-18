#ifndef TConstrainedLogLikelihood_H_seen
#define TConstrainedLikelihood_H_seen
// Generate a posterior for a likelihood that has constraints between
// different parametrers.  This example constrains the sum of the parameters
// and sets an estimated prior for each parameter value.
class TConstrainedLikelihood {
    std::vector<double> ExpectedValues;
    std::vector<double> PriorConstraints;
    double SummedValues;
    double SummedConstraint;
    
public:
    // Determine the number of dimensions.  This is where the dimensions are
    // defined, and everything else uses it.
    std::size_t GetDim() const {
        std::size_t dimensions = ExpectedValues.size();
        if (dimensions != PriorConstraints.size()) {
            std::cout << "IMPOSSIBLE DIMENSIONS Constrained Likelihood"
                      << std::endl;
            exit(1);
        }
        return dimensions;
    }

    // Calculate the log(likelihood).  The priors are set in the Init() method.
    double operator()(const Vector& point)  const {
        double logLikelihood = 0.0;


        // The sum should be constrained.
        double sum = 0.0;
        for (std::size_t i = 0; i<GetDim(); ++i) {
            sum += point[i];
        }
        sum = (sum-SummedValues)/SummedConstraint;
        logLikelihood -= 0.5*sum*sum;

        // The individual values should be constrained.
        for (std::size_t i = 0; i<GetDim(); ++i) {
            double v = point[i] - ExpectedValues[i];
            v /= PriorConstraints[i];
            logLikelihood -= 0.5*v*v;
        }
        
        return logLikelihood;
    }

    // Note that this needs to be the grad(log(Likelihood)).  It returns false
    // since the gradient is not calculated.
    bool operator() (Vector& g, const Vector& p) {
        return false;
    }

    // Set the priors for the likelihood.
    void Init() {
        // The sum of the parameters.
        SummedValues = 1902.0;
        // The uncertainty on the sum.
        SummedConstraint = 16.0;
        // 1
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 2
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 3
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 4
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 5
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 6
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 7
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 8
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 9
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 10
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 11
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 12
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 13
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 14
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 15
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 16
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 17
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 18
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 19
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 20
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 21
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 22
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 23
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 24
        ExpectedValues.push_back(76.0); PriorConstraints.push_back(76.0*0.08);
        // 25
        ExpectedValues.push_back(80.0); PriorConstraints.push_back(2.0);
    }
    
};
#endif
