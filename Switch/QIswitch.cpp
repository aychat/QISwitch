// tls.cpp : Defines the entry point for the console application.
//
#pragma warning( disable : 4996)

#include <stdio.h>
//#include <tchar.h>

#include "Modules/QIM_Vibronic.h"
//#include "Modules/QIM_Vibronic_ext.h"
#include "QI_Publisher.h"
#include "QI_Subscriber.h"
#include "QI_Liouvillian.h"
#include "operatorsInit.h"
#include "convert.h"

#include "propagation.inl"

using namespace util::convert; //defines $ function as universal converter of numbers to strings

void InitModel(QI_Model& m);
void InitParameter(QI_Model& m);

int main()
{
	QI_Model m;	//base object which collects all the information related to structure of Hilbert space and model parameters
	// ************** modules import
	QIM_Harmonic::Init(&m);	// Eigenenergies, creation, annihilation operators etc
	//QIM_Harmonic_EXT::Init(&m); // Franck-Condon transitions and their variations
	// ************** end of modules import

	InitModel(m);     // specification of the model geometry and declaration of model parameters
	InitParameter(m); // initialization of the parameters 
	// few extra propagation-related parameters
	m.Par() << "REAL max_prop_time";
	m.Par<double>("max_prop_time") = 150;
	m.Par() << "INT N_prop_iterations";
	m.Par<int>("N_prop_iterations") = 3000;
	m.Par() << "REAL pulse_duration";
	m.Par<double>("pulse_duration") = 10;
	m.Par() << "REAL pulse_strength";
	m.Par<double>("pulse_strength") = 10;
	m.Par() << "REAL max_energy";
	m.Par<double>("max_energy") = 32;

	//Specification of the system Hamiltonian
	QI_Matrix H_self(&m);  // Eigen (field-free) Hamiltonian
	QI_Matrix H_int(&m);   // Field-matter interaction part
	std::vector<QI_Matrix*> L_matrices;
	GetSelfHamiltonian(m, H_self);
	GetInteractionHamiltonian(m, H_int);
	GetRelaxationMatrices(m, L_matrices);

	//reindexation of the model by removing physically meaningless levels
	std::cout<< "Model size before indexing: "<<m.Space().Size()<<"\n";
	auto idx = m.Space().CreateIndex<QI_Topology::IndexType::RandomConstrained>();
	idx->Diagonalnit(H_self, -10000, m.Par<double>("max_energy"));
	std::cout<< "Model size after indexing: "<<m.Space().Size()<<"\n";
	//do not create observables before this point!

	// Definition of field-free and dressed parts of the Liouvillian
	QI_Liouvillian Lvn_0(&m), Lvn_i(&m);
	Lvn_0.AssignHamiltonian(H_self);
	std::for_each(L_matrices.begin(), L_matrices.end(), [&](QI_Matrix* L) {Lvn_0.AddLindbladOperator(*L); });
	Lvn_i.AssignHamiltonian(H_int);
	
	// Optimization of the numerical efficacy of propagator
	Lvn_0.Optimize(m.Space().Size(),m.Space().Size());
	Lvn_i.Optimize(m.Space().Size(), m.Space().Size());

	// Runge-Cutta propagation implementation
	double current_field_strength=0;
	int k = 0;
	double time = 0;
	double delta_t = 1.e-5;
	double precision = 1.e-8;
	auto advance = [&](QI_Observable* _output, QI_Observable* _input, double _prefactor, QI_Observable* _tmp)
	{
		Lvn_0.ApplyToState(_output, _input, _prefactor, _tmp);
		if(time<m.Par<double>("pulse_duration"))
			Lvn_i.ApplyToState(_output, _input, _prefactor*m.Par<double>("pulse_strength"), _tmp);
	};
	auto new_delta_t = [&](double _delta_t, double _inacc)
	{
		return delta_t*std::pow(1.1*precision / (_inacc + 0.1*precision), 0.2) / 1.05;
	};

	QI_Observable x(&m), x_old(&m),tmp(&m);
	x.SetZero();
	x(0, 0) = 1; //initial state;
	//x.SetIdentity();
	x_old = x;

	while ((++k < m.Par<int>("N_prop_iterations")) && (time < m.Par<double>("max_prop_time")))//for (int k = 0; k < 3000; k++)
	{
		auto inacc = MakeStep(x, delta_t,advance);
		delta_t = new_delta_t(delta_t, inacc);
		std::cout << delta_t << "*|" << time << "|"<<inacc<<"\n";
		while (inacc > precision)
		{
			delta_t = new_delta_t(delta_t, inacc);
			std::cout << delta_t << "-\n";
			x = x_old;
			inacc = MakeStep(x, delta_t, advance);
		};
		time += delta_t;
		x_old = x;
	}

	// clearning the memory
	std::for_each(L_matrices.begin(), L_matrices.end(), [](QI_Matrix* L) { delete L; });
    return 0;
};


void InitModel(QI_Model& m)
{
	//***************************************
	// specifying names for degrees of freedom (naming conventions are arbitrary, just avoid using ".", " " and ",")
	QI_DOF& d = m.DOF();
	d << "?e" << "?v0" << "?v1";
	
	//***************************************
	// specifying the structure of Configurational (Hilbert) space 
	//(the second argument is number of levels, 
	//the third is label for accessing the corresponding block, 
	//aka("label") is optional and should be used only to introduce a label for compound block)
	// + and * here imply direct product and direct sum
	m.Space() =  
		(d("?e", 1, "#E[e1]")*d("?v0", 1, "#E[e1_v0]")*d("?v1", 1, "#E[e1_v1]")).aka("#G")  //ground state
		+
		(d("?e", 1, "#E[e2]")*d("?v0", 10, "#E[e2_v0]")*d("?v1", 10, "#E[e2_v1]")).aka("#E") //excited state 
		+
		d("?e", 1, "#E[e3]")*d("?v0", 10, "#E[e3_v0]")*d("?v1", 10, "#E[e3_v1]") //final state
		+
		d("?e,?v0,?v1",1,"#D"); //drain state to emulate the population leakage out of the 3 electronic level subsystem

	//It is highly recommended to ensure that the total number of levels is chosen _e_v_e_n_. 
	//In the case of odd numbers program might
	//crush on some platforms and will always show ~1.5-2 times lower performance on the systems which allow for AVX

	//***************************************
	// Now we need to declare the model parameters. In is convenient to organize them into groups (publications)
	// according to their physical meaning. Then, specific subset of parameter can be accessed by "subscribing" to the 
	// appropriate groups. Next line defines the group names to be used.
	m.Pub() << "MOLECULAR_ENERGIES_AND_GEOMETRY"
		<< "ELECTRONIC_ENERGIES"
		<< "VIBRATIONAL_PROPERTIES"
		<< "ELECTRONIC_COUPLINGS"
		<< "VIBRATIONAL_COUPLINGS"
		<< "OPTICAL_EXCITATIONS"
		<< "QUANTUM_FRICTION"
		<< "DRAIN_DECAY";

	// Now we can introduce the model parameters as follows:
	m.Par().DefaultPub() <<= "MOLECULAR_ENERGIES_AND_GEOMETRY";
	m.Par().DefaultPub() << "ELECTRONIC_ENERGIES";
	m.Par()
		<< m("#E[e1]") //scope to which the following parameter applies, may be any logical expression composed of block labels introduced earlier (see below for examples)
		<< "std::multiplier			E[e1]"
		<< m("#E[e2]")
		<< "std::multiplier			E[e2]"
		<< m("#E[e3]")
		<< "std::multiplier			E[e3]";
	m.Par().DefaultPub() >> "ELECTRONIC_ENERGIES";

	m.Par().DefaultPub() << "VIBRATIONAL_PROPERTIES";
	m.Par() 
		<< m("?v0") 
			<< "harmonic::mass			mass0"
			<< "harmonic::frequency		w0"
			<< m("#E[e1_v0]")
			<< "harmonic::coord_eq		x_eq0[e1]"
			<< m("#E[e2_v0]")
			<< "harmonic::coord_eq		x_eq0[e2]"
			<< m("#E[e3_v0]")
			<< "harmonic::coord_eq		x_eq0[e3]"
		<< m("?v1")
			<< "harmonic::mass			mass1"
			<< "harmonic::frequency		w1"
			<< m("#E[e1_v1]")
			<< "harmonic::coord_eq		x_eq1[e1]"
			<< m("#E[e2_v1]")
			<< "harmonic::coord_eq		x_eq1[e2]"
			<< m("#E[e3_v1]")
			<< "harmonic::coord_eq		x_eq1[e3]";
	m.Par().DefaultPub() >> "VIBRATIONAL_PROPERTIES";

	m.Par().DefaultPub() <<= "VIBRATIONAL_COUPLINGS";
	m.Par()
		<< m("#E[e2]")
		<< "std::multiplier e2_v12"
		<< m("#E[e3]")
		<< "std::multiplier e3_v12";
	m.Par().DefaultPub() >> ("VIBRATIONAL_COUPLINGS");

	m.Par().DefaultPub() <<= "ELECTRONIC_COUPLINGS";
	m.Par()
		<< (m("#E[e2]") || m("#E[e3]")) // example of compound scope
		<< "std::multiplier e23";
	m.Par().DefaultPub() >> ("ELECTRONIC_COUPLINGS");

	m.Par().DefaultPub() <<= "OPTICAL_EXCITATIONS";
	m.Par()
		<< (m("#E[e1]") || m("#E[e2]")) // example of compound scope
		<< "std::multiplier field_strength";
	m.Par().DefaultPub() >> ("OPTICAL_EXCITATIONS");

	m.Par().DefaultPub() <<= "QUANTUM_FRICTION";
	m.Par() << (!m("?e")) << "harmonic::bondarian_friction_formula bondarian_velocity_dependence"; //applies for all vibrational degrees of freedom
	m.Par()
		<< (m("?v0") || m("?v1")) 
			<< "REAL friction_sign"
			<< "harmonic::formula friction_formula"
		<< m("?v0") << "std::multiplier	v0_friction"
		<< m("?v1") << "std::multiplier	v1_friction";
	m.Par().DefaultPub() >> "QUANTUM_FRICTION";

	m.Par().DefaultPub() <<= "DRAIN_DECAY";
		m.Par() << (m("#D") || !m("#D")) // i.e apply it to everything
			<< "std::multiplier	drain_decay";
	m.Par().DefaultPub() >> "DRAIN_DECAY";
};

void InitParameter(QI_Model& m)
{
	// here we assign actual values to all parameters;

	auto pc = [&](std::string s)->std::complex<double>&
	{
		return m.Par<std::complex<double> >(s); // generic syntax for parameter access 
	};
	auto pd = [&](std::string s)->double&
	{
		return m.Par<double>(s); 
	};

	m.Par<double>("std_constant.hbar") = 1.; //!!! rescaling the units here

	pc("E[e1]") = 0.001; // equivalent to m.Par<std::complex<double> >("E[el]") = 0.;
	pc("E[e2]") = 20;
	pc("E[e3]") = 14;

	pd("mass0") = 1;
	pd("w0") = 1;
	pd("x_eq0[e1]") = 0.001;
	pd("x_eq0[e2]") = 1;
	pd("x_eq0[e3]") = 2;

	pd("mass1") = 1;
	pd("w1") = 1;
	pd("x_eq1[e1]") = 0.001;
	pd("x_eq1[e2]") = 1;
	pd("x_eq1[e3]") = 2;

	pc("e2_v12") = 0.5;
	pc("e3_v12") = -0.5;

	pc("e23") = 1;

	pc("field_strength") = 1.;

	m.Par<std::string>("bondarian_velocity_dependence") = "rect(kappa_bdr*hbar/2,20.0,p_bdr)*(((p_bdr-kappa_bdr*hbar/2)*(p_bdr-kappa_bdr*hbar/2))^0.25)";
	m.Par<std::string>("friction_formula")= "Bondarian(2,par(friction_sign))";
	pc("v0_friction") = 0.001;
	pc("v1_friction") = 0.001;

	pc("drain_decay") = 0.0005;
};

