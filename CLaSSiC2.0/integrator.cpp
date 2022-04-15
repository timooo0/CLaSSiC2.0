#include "integrator.hpp"
#include <iostream>
#include <cmath>
#include <immintrin.h>

void Integrator::calculateEffectiveField(std::vector<std::vector<int>> &neighbours, Vector3D &spin)
{
	/*
	Calculates the effective by combining: external field, anisotropy and nearest neighbour interaction
	*/
	float jx;
	float jy;
	float jz;

/*	
	__m256 mFieldX;
	__m256 mFieldY;
	__m256 mFieldZ;
	__m256 mSpinZ;
	__m256 mAnisAxis = _mm256_set1_ps(constants::anisotropyAxis);
	__m256 mAnisPlane = _mm256_set1_ps(constants::anisotropyPlane);

	for (int i = 0; i < constants::nAtoms; i+=8)
	{
		// External field
		// memcpy(&effectiveField.x[3 * i], &constants::magneticField, sizeoff(constants::magneticField)); // set to external field

		mFieldX = _mm256_loadu_ps(&effectiveField.x[i]);
		mFieldY = _mm256_loadu_ps(&effectiveField.y[i]);
		mFieldZ = _mm256_set1_ps(constants::magneticField[2]);
		mSpinZ = _mm256_loadu_ps(&spin.z[i]);

		//effectiveField.z[i] = constants::magneticField[2];

		// Anisotropy
		mFieldZ = _mm256_fmadd_ps(mAnisAxis, mSpinZ, mFieldZ);
		mFieldZ = _mm256_fmsub_ps(mAnisPlane, mSpinZ, mFieldZ);
		// effectiveField.z[i] += constants::anisotropyAxis * spin.z[i];
		// effectiveField.z[i] -= constants::anisotropyPlane * spin.z[i];

		_mm256_storeu_ps(&effectiveField.z[i], _mm256_fmadd_ps(mSpinZ, _mm256_sub_ps(mAnisAxis, mAnisPlane), _mm256_set1_ps(constants::magneticField[2])));
		//Stabilizer for kagome lattice by adding anisotropy (1 T) alligned with GS
		// if (constants::geometry==4 && constants::stabilize){
		// 	effectiveField.x[i] += stabilizerField.x[i];
		// 	effectiveField.y[i] += stabilizerField.y[i];
		// }
	}*/
	
	for (int i = 0; i < constants::nAtoms; i++)
	{
		// External field
		// memcpy(&effectiveField.x[3 * i], &constants::magneticField, sizeoff(constants::magneticField)); // set to external field

		effectiveField.z[i] = constants::magneticField[2];

		// Anisotropy
		effectiveField.z[i] += constants::anisotropyAxis * spin.z[i];
		effectiveField.z[i] -= constants::anisotropyPlane * spin.z[i];

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
	/*
	__m256 mSpinX;
	__m256 mSpinY;
	__m256 mSpinZ;
	__m256 mFieldX;
	__m256 mFieldY;
	__m256 mFieldZ;
	__m256 mGammaDt = _mm256_set1_ps(-constants::gamma * constants::dt);
	__m256 mLambda = _mm256_set1_ps(-constants::lambda);

	for (int i = 0; i < constants::nAtoms; i+=8)
	{
		mSpinX = _mm256_loadu_ps(&spin.x[i]);
		mSpinY = _mm256_loadu_ps(&spin.y[i]);
		mSpinZ = _mm256_loadu_ps(&spin.z[i]);
		mFieldX = _mm256_loadu_ps(&effectiveField.x[i]);
		mFieldY = _mm256_loadu_ps(&effectiveField.y[i]);
		mFieldZ = _mm256_loadu_ps(&effectiveField.z[i]);

		_mm256_storeu_ps(&rk.x[i], 
		 _mm256_mul_ps(mGammaDt, 
		 _mm256_add_ps(_mm256_sub_ps(_mm256_mul_ps(mSpinY, mFieldZ), _mm256_mul_ps(mSpinZ, mFieldY)), 
		   _mm256_mul_ps(mLambda,
		    _mm256_add_ps(
		     _mm256_sub_ps(_mm256_mul_ps(_mm256_mul_ps(mFieldY, mSpinX), mSpinY), _mm256_mul_ps(_mm256_mul_ps(mFieldX, mSpinY), mSpinY)), 
		     _mm256_sub_ps( _mm256_mul_ps(_mm256_mul_ps(mFieldZ, mSpinX), mSpinZ), _mm256_mul_ps(_mm256_mul_ps(mFieldX, mSpinZ), mSpinZ))
		    )
		   )
		  )
		 )
		);
		
		_mm256_storeu_ps(&rk.y[i], 
		 _mm256_mul_ps(mGammaDt, 
		  _mm256_add_ps(_mm256_sub_ps(_mm256_mul_ps(mSpinZ, mFieldX), _mm256_mul_ps(mSpinX, mFieldZ)), 
		   _mm256_mul_ps(mLambda,
		    _mm256_add_ps(
		     _mm256_sub_ps(_mm256_mul_ps(_mm256_mul_ps(mFieldX, mSpinX), mSpinY), _mm256_mul_ps(_mm256_mul_ps(mFieldY, mSpinX), mSpinX)),
		     _mm256_sub_ps( _mm256_mul_ps(_mm256_mul_ps(mFieldZ, mSpinY), mSpinZ), _mm256_mul_ps(_mm256_mul_ps(mFieldY, mSpinZ), mSpinZ))
		    )
		   )
		  )
		 )
		);

		_mm256_storeu_ps(&rk.z[i],
		 _mm256_mul_ps(mGammaDt,
		  _mm256_add_ps(_mm256_sub_ps(_mm256_mul_ps(mSpinX, mFieldY), _mm256_mul_ps(mSpinY, mFieldX)), 
		   _mm256_mul_ps(mLambda,
		    _mm256_add_ps(
		     _mm256_sub_ps(_mm256_mul_ps(_mm256_mul_ps(mFieldY, mSpinY), mSpinZ), _mm256_mul_ps(_mm256_mul_ps(mFieldZ, mSpinY), mSpinY)),
		    _mm256_sub_ps( _mm256_mul_ps(_mm256_mul_ps(mFieldX, mSpinX), mSpinZ), _mm256_mul_ps(_mm256_mul_ps(mFieldZ, mSpinX), mSpinX))
			)
		   )
		  )
		 )
		);
	} */
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
	}
}

void Integrator::integrate(std::vector<std::vector<int>> &neighbours, Vector3D &spin, std::vector<float> &randomField)
{
	//More hardcore? Does Euler work?
	evaluate(neighbours, spin);

	//Update rkPos
	/*
	__m256 a;
	__m256 b;
	__m256 halfVec = _mm256_set1_ps(0.5f);
	__m256 result;
	for (int i = 0; i < constants::nAtoms; i+=8)
	{
		a = _mm256_loadu_ps(&spin.x[i]);
		b = _mm256_loadu_ps(&rk.x[i]);
		result = _mm256_fmadd_ps(a, halfVec, b);
		_mm256_storeu_ps(&rkPos.x[i], result);
	}
	for (int i = 0; i < constants::nAtoms; i+=8)
	{
		a = _mm256_loadu_ps(&spin.y[i]);
		b = _mm256_loadu_ps(&rk.y[i]);
		result = _mm256_fmadd_ps(a, halfVec, b);
		_mm256_storeu_ps(&rkPos.y[i], result);
	}
	for (int i = 0; i < constants::nAtoms; i+=8)
	{
		a = _mm256_loadu_ps(&spin.z[i]);
		b = _mm256_loadu_ps(&rk.z[i]);
		result = _mm256_fmadd_ps(a, halfVec, b);
		_mm256_storeu_ps(&rkPos.z[i], result);
	}*/
	for (int i = 0; i < constants::nAtoms; i++)
	{
		rkPos.x[i] = spin.x[i] + rk.x[i] / 2.;
		rkPos.y[i] = spin.y[i] + rk.y[i] / 2.;
		rkPos.z[i] = spin.z[i] + rk.z[i] / 2.;
	}
	evaluate(neighbours, spin);

/*
	__m256 mRandomX;
	__m256 mRandomY;
	__m256 mRandomZ;
	__m256 mSpinX;
	__m256 mSpinY;
	__m256 mSpinZ;
	__m256 mRkX;
	__m256 mRkY;
	__m256 mRkZ;
	__m256 mGamma = _mm256_set1_ps(constants::gamma);

	for (int i = 0; i < constants::nAtoms; i+=8)
	{
		mSpinX = _mm256_loadu_ps(&spin.x[i]);
		mSpinY = _mm256_loadu_ps(&spin.y[i]);
		mSpinZ = _mm256_loadu_ps(&spin.z[i]);
		mRkX = _mm256_loadu_ps(&rk.x[i]);
		mRkY = _mm256_loadu_ps(&rk.y[i]);
		mRkZ = _mm256_loadu_ps(&rk.z[i]);
		mRandomX = _mm256_loadu_ps(&randomField[3*i]);
		mRandomY = _mm256_loadu_ps(&randomField[3*i+8]);
		mRandomZ = _mm256_loadu_ps(&randomField[3*i+16]);

		_mm256_storeu_ps(&spin.x[i], _mm256_add_ps(mSpinX, _mm256_add_ps(mRkX, _mm256_mul_ps(mGamma, _mm256_sub_ps(_mm256_mul_ps(mSpinY, mRandomZ), _mm256_mul_ps(mSpinZ, mRandomY))))));
		_mm256_storeu_ps(&spin.y[i], _mm256_add_ps(mSpinY, _mm256_add_ps(mRkY, _mm256_mul_ps(mGamma, _mm256_sub_ps(_mm256_mul_ps(mSpinZ, mRandomX), _mm256_mul_ps(mSpinX, mRandomZ))))));
		_mm256_storeu_ps(&spin.z[i], _mm256_add_ps(mSpinZ, _mm256_add_ps(mRkZ, _mm256_mul_ps(mGamma, _mm256_sub_ps(_mm256_mul_ps(mSpinX, mRandomY), _mm256_mul_ps(mSpinY, mRandomX))))));
	}*/
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
}

float Integrator::calculateEnergy(Vector3D &spin){
	float totalEnergy = 0;
	for (int i=0;i<3*constants::nAtoms;i++){
		totalEnergy += spin.x[i] * effectiveField.x[i];
		totalEnergy += spin.y[i] * effectiveField.y[i];
		totalEnergy += spin.z[i] * effectiveField.z[i]; 
	}
	return -constants::gFactor*constants::bohrMagneton*totalEnergy;
}
float Integrator::dotProduct(const float a[3], float b[3])
{
	float result = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	return result;
}

Vector3D* Integrator::getEffectiveField(){
	return &effectiveField;
};
void Integrator::setStabilizerField(){
	float iterAngle;
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
			for (float i = 0; i < constants::nAtoms; i++){
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

[[nodiscard]] inline static __m256 crossElement( __m256 const& a, __m256 const& b, __m256 const& c, __m256 const& d ) {
    return _mm256_fmadd_ps(a, b, _mm256_add_ps(c, d));
	}