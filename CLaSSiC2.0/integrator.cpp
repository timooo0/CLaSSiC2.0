#include "integrator.hpp"
#include <iostream>
#include <cmath>
#include <immintrin.h>

void Integrator::calculateEffectiveField(std::vector<std::vector<int>> &neighbours, Vector3D &spin)
{
	/*
	Calculates the effective by combining: external field, anisotropy and nearest neighbour interaction
	*/
	double jx;
	double jy;
	double jz;

	
	__m256d mFieldX;
	__m256d mFieldY;
	__m256d mFieldZ;
	__m256d mSpinZ;
	__m256d mAnisAxis = _mm256_set1_pd(constants::anisotropyAxis);
	__m256d mAnisPlane = _mm256_set1_pd(constants::anisotropyPlane);

	for (int i = 0; i < constants::nAtoms; i+=4)
	{
		// External field
		// memcpy(&effectiveField.x[3 * i], &constants::magneticField, sizeoff(constants::magneticField)); // set to external field

		mFieldX = _mm256_loadu_pd(&effectiveField.x[i]);
		mFieldY = _mm256_loadu_pd(&effectiveField.y[i]);
		mFieldZ = _mm256_set1_pd(constants::magneticField[2]);
		mSpinZ = _mm256_loadu_pd(&spin.z[i]);

		//effectiveField.z[i] = constants::magneticField[2];

		// Anisotropy
		mFieldZ = _mm256_fmadd_pd(mAnisAxis, mSpinZ, mFieldZ);
		mFieldZ = _mm256_fmsub_pd(mAnisPlane, mSpinZ, mFieldZ);
		// effectiveField.z[i] += constants::anisotropyAxis * spin.z[i];
		// effectiveField.z[i] -= constants::anisotropyPlane * spin.z[i];

		_mm256_storeu_pd(&effectiveField.z[i], _mm256_fmadd_pd(mSpinZ, _mm256_sub_pd(mAnisAxis, mAnisPlane), _mm256_set1_pd(constants::magneticField[2])));
		//Stabilizer for kagome lattice by adding anisotropy (1 T) alligned with GS
		// if (constants::geometry==4 && constants::stabilize){
		// 	effectiveField.x[i] += stabilizerField.x[i];
		// 	effectiveField.y[i] += stabilizerField.y[i];
		// }
	}
	
	for (int i = 0; i < constants::nAtoms; i++)
	{
		// External field
		// memcpy(&effectiveField.x[3 * i], &constants::magneticField, sizeoff(constants::magneticField)); // set to external field
		/*
		effectiveField.z[i] = constants::magneticField[2];

		// Anisotropy
		effectiveField.z[i] += constants::anisotropyAxis * spin.z[i];
		effectiveField.z[i] -= constants::anisotropyPlane * spin.z[i];
		*/
		// Exchange
		jx = 0;
		jy = 0;
		jz = 0;
		// Nearest neighbours
		for (int j : neighbours[i])
		{
			jx += spin.x[j];	
			jy += spin.y[j];
			jz += spin.z[j];
		}
		effectiveField.x[i] += constants::exchangePrefactor * jx; // x
		effectiveField.y[i] += constants::exchangePrefactor * jy; // y
		effectiveField.z[i] += constants::exchangePrefactor * jz; // z
	}
		//Stabilizer for kagome lattice by adding anisotropy (1 T) alligned with GS
		// if (constants::geometry==4 && constants::stabilize){
		// 	effectiveField.x[i] += stabilizerField.x[i];
		// 	effectiveField.y[i] += stabilizerField.y[i];
		// }

}

void Integrator::evaluate(std::vector<std::vector<int>> &neighbours, Vector3D &spin)
{
	/*
	function evaluation
	*/
	calculateEffectiveField(neighbours, spin);
	
	__m256d mSpinX;
	__m256d mSpinY;
	__m256d mSpinZ;
	__m256d mFieldX;
	__m256d mFieldY;
	__m256d mFieldZ;
	__m256d mGammaDt = _mm256_set1_pd(-constants::gamma * constants::dt);
	__m256d mLambda = _mm256_set1_pd(-constants::lambda);

	for (int i = 0; i < constants::nAtoms; i+=4)
	{
		mSpinX = _mm256_loadu_pd(&spin.x[i]);
		mSpinY = _mm256_loadu_pd(&spin.y[i]);
		mSpinZ = _mm256_loadu_pd(&spin.z[i]);
		mFieldX = _mm256_loadu_pd(&effectiveField.x[i]);
		mFieldY = _mm256_loadu_pd(&effectiveField.y[i]);
		mFieldZ = _mm256_loadu_pd(&effectiveField.z[i]);

		_mm256_storeu_pd(&rk.x[i], 
		 _mm256_mul_pd(mGammaDt, 
		 _mm256_add_pd(_mm256_sub_pd(_mm256_mul_pd(mSpinY, mFieldZ), _mm256_mul_pd(mSpinZ, mFieldY)), 
		   _mm256_mul_pd(mLambda,
		    _mm256_add_pd(
		     _mm256_sub_pd(_mm256_mul_pd(_mm256_mul_pd(mFieldY, mSpinX), mSpinY), _mm256_mul_pd(_mm256_mul_pd(mFieldX, mSpinY), mSpinY)), 
		     _mm256_sub_pd( _mm256_mul_pd(_mm256_mul_pd(mFieldZ, mSpinX), mSpinZ), _mm256_mul_pd(_mm256_mul_pd(mFieldX, mSpinZ), mSpinZ))
		    )
		   )
		  )
		 )
		);
		
		_mm256_storeu_pd(&rk.y[i], 
		 _mm256_mul_pd(mGammaDt, 
		  _mm256_add_pd(_mm256_sub_pd(_mm256_mul_pd(mSpinZ, mFieldX), _mm256_mul_pd(mSpinX, mFieldZ)), 
		   _mm256_mul_pd(mLambda,
		    _mm256_add_pd(
		     _mm256_sub_pd(_mm256_mul_pd(_mm256_mul_pd(mFieldX, mSpinX), mSpinY), _mm256_mul_pd(_mm256_mul_pd(mFieldY, mSpinX), mSpinX)),
		     _mm256_sub_pd( _mm256_mul_pd(_mm256_mul_pd(mFieldZ, mSpinY), mSpinZ), _mm256_mul_pd(_mm256_mul_pd(mFieldY, mSpinZ), mSpinZ))
		    )
		   )
		  )
		 )
		);

		_mm256_storeu_pd(&rk.z[i],
		 _mm256_mul_pd(mGammaDt,
		  _mm256_add_pd(_mm256_sub_pd(_mm256_mul_pd(mSpinX, mFieldY), _mm256_mul_pd(mSpinY, mFieldX)), 
		   _mm256_mul_pd(mLambda,
		    _mm256_add_pd(
		     _mm256_sub_pd(_mm256_mul_pd(_mm256_mul_pd(mFieldY, mSpinY), mSpinZ), _mm256_mul_pd(_mm256_mul_pd(mFieldZ, mSpinY), mSpinY)),
		    _mm256_sub_pd( _mm256_mul_pd(_mm256_mul_pd(mFieldX, mSpinX), mSpinZ), _mm256_mul_pd(_mm256_mul_pd(mFieldZ, mSpinX), mSpinX))
			)
		   )
		  )
		 )
		);
	}
	/*
	for (int i = 0; i < constants::nAtoms; i++)
	{
		rk.x[i] = -constants::gamma * constants::dt * (spin.y[i] * effectiveField.z[i] - spin.z[i] * effectiveField.y[i]);
		// - constants::lambda * 
		//  (effectiveField.y[i] * spin.x[i] * spin.y[i]
		//   - effectiveField.x[i] * spin.y[i] * spin.y[i]
		//   + effectiveField.z[i] * spin.x[i] * spin.z[i]
		//    - effectiveField.x[i] * spin.z[i] * spin.z[i]));
		
		rk.y[i] = -constants::gamma * constants::dt * (spin.z[i] * effectiveField.x[i] - spin.x[i] * effectiveField.z[i]);
		// - constants::lambda * 
		//  (-effectiveField.y[i] * spin.x[i] * spin.x[i]
		//   + effectiveField.x[i] * spin.x[i] * spin.y[i]
		//   + effectiveField.z[i] * spin.y[i] * spin.z[i]
		//    - effectiveField.y[i] * spin.z[i] * spin.z[i]));


		rk.z[i] = -constants::gamma * constants::dt * (spin.x[i] * effectiveField.y[i] - spin.y[i] * effectiveField.x[i]); 
		// - constants::lambda * 
		// (-effectiveField.z[i] * spin.x[i] * spin.x[i]
		//  - effectiveField.z[i] * spin.y[i] * spin.y[i]
		//  + effectiveField.x[i] * spin.x[i] * spin.z[i] 
		//  + effectiveField.y[i] * spin.y[i] * spin.z[i]));
	}*/
}

void Integrator::integrate(std::vector<std::vector<int>> &neighbours, Vector3D &spin, std::vector<double> &randomField)
{
	//More hardcore? Does Euler work?
	evaluate(neighbours, spin);
	//Update rkPos
	
	__m256d a;
	__m256d b;
	__m256d halfVec = _mm256_set1_pd(0.5f);
	__m256d result;
	for (int i = 0; i < constants::nAtoms; i+=4)
	{
		a = _mm256_loadu_pd(&rk.x[i]);
		b = _mm256_loadu_pd(&spin.x[i]);
		result = _mm256_fmadd_pd(a, halfVec, b);
		_mm256_storeu_pd(&rkPos.x[i], result);
	}
	for (int i = 0; i < constants::nAtoms; i+=4)
	{
		a = _mm256_loadu_pd(&rk.y[i]);
		b = _mm256_loadu_pd(&spin.y[i]);
		result = _mm256_fmadd_pd(a, halfVec, b);
		_mm256_storeu_pd(&rkPos.y[i], result);
	}
	for (int i = 0; i < constants::nAtoms; i+=4)
	{
		a = _mm256_loadu_pd(&rk.z[i]);
		b = _mm256_loadu_pd(&spin.z[i]);
		result = _mm256_fmadd_pd(a, halfVec, b);
		_mm256_storeu_pd(&rkPos.z[i], result);
	}
	/*
	for (int i = 0; i < constants::nAtoms; i++)
	{
		rkPos.x[i] = spin.x[i] + rk.x[i] / 2.;
		rkPos.y[i] = spin.y[i] + rk.y[i] / 2.;
		rkPos.z[i] = spin.z[i] + rk.z[i] / 2.;
	}
	*/
	evaluate(neighbours, rkPos);


	__m256d mRandomX;
	__m256d mRandomY;
	__m256d mRandomZ;
	__m256d mSpinX;
	__m256d mSpinY;
	__m256d mSpinZ;
	__m256d mRkX;
	__m256d mRkY;
	__m256d mRkZ;
	__m256d mGamma = _mm256_set1_pd(constants::gamma);

	for (int i = 0; i < constants::nAtoms; i+=4)
	{
		mSpinX = _mm256_loadu_pd(&spin.x[i]);
		mSpinY = _mm256_loadu_pd(&spin.y[i]);
		mSpinZ = _mm256_loadu_pd(&spin.z[i]);
		mRkX = _mm256_loadu_pd(&rk.x[i]);
		mRkY = _mm256_loadu_pd(&rk.y[i]);
		mRkZ = _mm256_loadu_pd(&rk.z[i]);
		mRandomX = _mm256_loadu_pd(&randomField[3*i]);
		mRandomY = _mm256_loadu_pd(&randomField[3*i+4]);
		mRandomZ = _mm256_loadu_pd(&randomField[3*i+8]);

		_mm256_storeu_pd(&spin.x[i], _mm256_add_pd(mSpinX, _mm256_add_pd(mRkX, _mm256_mul_pd(mGamma, _mm256_sub_pd(_mm256_mul_pd(mSpinY, mRandomZ), _mm256_mul_pd(mSpinZ, mRandomY))))));
		_mm256_storeu_pd(&spin.y[i], _mm256_add_pd(mSpinY, _mm256_add_pd(mRkY, _mm256_mul_pd(mGamma, _mm256_sub_pd(_mm256_mul_pd(mSpinZ, mRandomX), _mm256_mul_pd(mSpinX, mRandomZ))))));
		_mm256_storeu_pd(&spin.z[i], _mm256_add_pd(mSpinZ, _mm256_add_pd(mRkZ, _mm256_mul_pd(mGamma, _mm256_sub_pd(_mm256_mul_pd(mSpinX, mRandomY), _mm256_mul_pd(mSpinY, mRandomX))))));
	}
	/*
	double temperatureSpin[3];
	for (int i = 0; i < constants::nAtoms; i++)
	{
		temperatureSpin[0] = constants::gamma * (spin.y[i] * randomField[3 * i + 2] - spin.z[i] * randomField[3 * i + 1]);
		temperatureSpin[1] = constants::gamma * (spin.z[i] * randomField[3 * i] - spin.x[i] * randomField[3 * i + 2]);
		temperatureSpin[2] = constants::gamma * (spin.x[i] * randomField[3 * i + 1] - spin.y[i] * randomField[3 * i]);

		spin.x[i] += rk.x[i];// + temperatureSpin[0];
		spin.y[i] += rk.y[i];// + temperatureSpin[1];
		spin.z[i] += rk.z[i];// + temperatureSpin[2];
	}
	*/
}

double Integrator::calculateEnergy(Vector3D &spin){
	double totalEnergy = 0;
	for (int i=0;i<3*constants::nAtoms;i++){
		totalEnergy += spin.x[i] * effectiveField.x[i];
		totalEnergy += spin.y[i] * effectiveField.y[i];
		totalEnergy += spin.z[i] * effectiveField.z[i]; 
	}
	return -constants::gFactor*constants::bohrMagneton*totalEnergy;
}
double Integrator::dotProduct(const double a[3], double b[3])
{
	double result = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	return result;
}

Vector3D* Integrator::getEffectiveField(){
	return &effectiveField;
};
void Integrator::setStabilizerField(){
	double iterAngle;
	switch (constants::spinInit) {
		case 7: 
			for (int i = 0; i < constants::nAtoms; i++){
				iterAngle = 1./2.*constants::pi-2./3.*constants::pi*(i/3)+2./3.*constants::pi*i;
				stabilizerField.x[i] = std::cos(iterAngle);
				stabilizerField.y[i] = std::sin(iterAngle);
				stabilizerField.z[i] = 0;
			}
			break;
		case 8:
			for (double i = 0; i < constants::nAtoms; i++){
			iterAngle = constants::pi*7./6.+2./3.*constants::pi*i;
			stabilizerField.x[i] = std::cos(iterAngle);
			stabilizerField.y[i] = std::sin(iterAngle);
			stabilizerField.z[i] = 0;
		}
	}
	for (int i=0;i<9;i++){
		std::cout << stabilizerField.x[i] << ", " << stabilizerField.y[i] << std::endl;
	}

}

[[nodiscard]] inline static __m256d crossElement( __m256d const& a, __m256d const& b, __m256d const& c, __m256d const& d ) {
    return _mm256_fmadd_pd(a, b, _mm256_add_pd(c, d));
	}