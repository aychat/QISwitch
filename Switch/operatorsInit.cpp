#pragma warning( disable : 4996)

#include "operatorsInit.h"
#include <string>
#include <ios>
#include "convert.h"

#include "QI_Publisher.h"
#include "QI_Subscriber.h"
#include "QI_Matrix.h"

using namespace util::convert;

void GetSelfHamiltonian(QI_Model& m, QI_Matrix& H_self)
{
	QI_MatrixSchema H_schema =
		m["#E[e1],#E[e1]"] * m["#E[e1_v0],#E[e1_v0]"] * m["#E[e1_v1],#E[e1_v1]"]
		+ m["#E[e2],#E[e2]"] * m["#E[e2_v0],#E[e2_v0]"] * m["#E[e2_v1],#E[e2_v1]"]
		+ m["#E[e3],#E[e3]"] * m["#E[e3_v0],#E[e3_v0]"] * m["#E[e3_v1],#E[e3_v1]"];
	
	// electronic energies
	{
		QI_Matrix H_tmp(&m);
		H_tmp.Subscriptions() = (m.Pub()("MOLECULAR_ENERGIES_AND_GEOMETRY") || m.Pub()("default") || m.Pub()("temp"));
		H_tmp.Schema() = H_schema;
		H_tmp.Alg(!m("?e")) = $$("std::IdentityMatrix");
		H_tmp.Alg(m("?e")) = $$("std::ScaledIdentityMatrix");
		H_tmp.Init();
		H_self = H_tmp;
	}

	// vibrational energies
	{
		QI_Matrix H_tmp(&m);
		H_tmp.Subscriptions() = (m.Pub()("MOLECULAR_ENERGIES_AND_GEOMETRY") || m.Pub()("default") || m.Pub()("temp"));
		H_tmp.Schema() = H_schema;
		H_tmp.Alg(!m("?v0")) = $$("std::IdentityMatrix");
		H_tmp.Alg(m("?v0")) = $$("harmonic::auto");
		H_tmp.Init();
		H_self += H_tmp;
	}
	{
		QI_Matrix H_tmp(&m);
		H_tmp.Subscriptions() = (m.Pub()("MOLECULAR_ENERGIES_AND_GEOMETRY") || m.Pub()("default") || m.Pub()("temp"));
		H_tmp.Schema() = H_schema;
		H_tmp.Alg(!m("?v1")) = $$("std::IdentityMatrix");
		H_tmp.Alg(m("?v1")) = $$("harmonic::auto");
		H_tmp.Init();
		H_self += H_tmp;
	}

	H_schema =
		m["#E[e2],#E[e2]"] * m["#E[e2_v0],#E[e2_v0]"] * m["#E[e2_v1],#E[e2_v1]"]
		+ m["#E[e3],#E[e3]"] * m["#E[e3_v0],#E[e3_v0]"] * m["#E[e3_v1],#E[e3_v1]"];
	// vibrational couplings
	{
		QI_Matrix H_tmp(&m, QI_Matrix::InitFlags::Hermitian);
		H_tmp.Subscriptions() = m.Pub()("VIBRATIONAL_COUPLINGS") || m.Pub()("VIBRATIONAL_PROPERTIES") || m.Pub()("default") || m.Pub()("temp");
		H_tmp.Alg(m("?e")) = $$("std::ScaledIdentityMatrix");
		H_tmp.Alg(m("?v0") || m("?v1")) = $$("harmonic::position_rel");
		H_tmp.Schema() = H_schema;
		H_tmp.Init();
		H_self += H_tmp;
	};
	// electronic couplings
	H_schema =
		m["#E[e2],#E[e3]"] * m["#E[e2_v0],#E[e3_v0]"] * m["#E[e2_v1],#E[e3_v1]"];
	{
		QI_Matrix H_tmp(&m, QI_Matrix::InitFlags::Hermitian);
		H_tmp.Subscriptions() = m.Pub()("ELECTRONIC_COUPLINGS") || m.Pub()("VIBRATIONAL_PROPERTIES") || m.Pub()("default") || m.Pub()("temp");
		H_tmp.Alg(m("?e")) = $$("std::ConstantMatrix");
		H_tmp.Alg(m("?v0") || m("?v1")) = $$("harmonic::auto");
		H_tmp.Schema() = H_schema;
		H_tmp.Init();
		H_self += H_tmp;
	};
};

void GetInteractionHamiltonian(QI_Model& m, QI_Matrix& H_int)
{
	// electronic couplings
	QI_MatrixSchema H_schema =
		m["#E[e1],#E[e2]"] * m["#E[e1_v0],#E[e2_v0]"] * m["#E[e1_v1],#E[e2_v1]"];
	{
		QI_Matrix H_tmp(&m, QI_Matrix::InitFlags::Hermitian);
		H_tmp.Subscriptions() = m.Pub()("OPTICAL_EXCITATIONS") || m.Pub()("VIBRATIONAL_PROPERTIES") || m.Pub()("default") || m.Pub()("temp");
		H_tmp.Alg(m("?e")) = $$("std::ScaledIdentityMatrix");//("std::ConstantMatrix");
		H_tmp.Alg(m("?v0") || m("?v1")) = $$("harmonic::auto");
		H_tmp.Schema() = H_schema;
		H_tmp.Init();
		H_int = H_tmp;
	};
};

void GetRelaxationMatrices(QI_Model& m, std::vector<QI_Matrix*>& L_matrices)
{
	//quantum friction
	QI_MatrixSchema L_schema =
		m["#E[e1],#E[e1]"] * m["#E[e1_v0],#E[e1_v0]"] * m["#E[e1_v1],#E[e1_v1]"]
		+ m["#E[e2],#E[e2]"] * m["#E[e2_v0],#E[e2_v0]"] * m["#E[e2_v1],#E[e2_v1]"]
		+ m["#E[e3],#E[e3]"] * m["#E[e3_v0],#E[e3_v0]"] * m["#E[e3_v1],#E[e3_v1]"];
	m.Par<double>("friction_sign") = 1;
	while (m.Par<double>("friction_sign") > -2)
	{
		{
			QI_Matrix* L_tmp = new QI_Matrix(&m);
	 		L_tmp->Schema() = L_schema;
			L_tmp->Subscriptions() = m.Pub()("QUANTUM_FRICTION") || m.Pub()("VIBRATIONAL_PROPERTIES") || m.Pub()("default") || m.Pub()("temp");
			L_tmp->Alg(!m("?v0")) = $$("std::IdentityMatrix");
			L_tmp->Alg(m("?v0")) = $$("harmonic::fromFormula");
			L_tmp->Init();
			L_matrices.push_back(L_tmp);
		}
		{
			QI_Matrix* L_tmp = new QI_Matrix(&m);
			L_tmp->Schema() = L_schema;
			L_tmp->Subscriptions() = m.Pub()("QUANTUM_FRICTION") || m.Pub()("VIBRATIONAL_PROPERTIES") || m.Pub()("default") || m.Pub()("temp");
			L_tmp->Alg(!m("?v1")) = $$("std::IdentityMatrix");
			L_tmp->Alg(m("?v1")) = $$("harmonic::fromFormula");
			L_tmp->Init();
			L_matrices.push_back(L_tmp);
		}
		m.Par<double>("friction_sign") -= 2;
	}
	
	//drain relaxation
	{	
		auto Add_Drain = [&](std::string drain, std::string state)
		{
			auto l_drain1_schema = m[drain + "," + state];
			auto l_drain1_subscriptions = m.Pub()("DRAIN_DECAY") || m.Pub()("default") || m.Pub()("temp");
			QI::local_idx dim = m.Space(state).Size();
			MatrixXc_r_hsp rel_block(1, dim);
			for (int i = 0; i < dim; i++)
			{
				rel_block.setZero();
				rel_block.coeffRef(0, i) = 1.;

				QI_Matrix* L_tmp = new QI_Matrix(&m);
				L_tmp->Schema() = L_schema;
				L_tmp->Schema() = l_drain1_schema; 
				L_tmp->Subscriptions() = l_drain1_subscriptions;
				L_tmp->Alg(m(state) || m(drain)) = $$("std::CustomSparseMatrix");
				L_tmp->Par("std::hypersparse_matrix", m(state) || m(drain)) = rel_block;
				L_tmp->Init();
				L_matrices.push_back(L_tmp);
			}
		};
		Add_Drain("#D", "#E");
	};
}