//
// Created by zhouxuguang on 2026/10/17.
// 快速傅里叶变换（FFT）
//

#include "FFT.h"

uint32_t ReverseBitsLoop(uint32_t x, uint32_t bits)
{
	if (bits == 0)
		return 0;

	uint32_t result = 0;
	for (uint32_t i = 0; i < bits; ++i)
	{
		result = (result << 1) | (x & 1u);
		x >>= 1;
	}
	return result;
}

FFT::FFT(int n)
{
    mN = n;
    mTw.resize(mN / 2);
	buildTwiddle();

	// 计算蝶形运算的级数
	int stages = 0;
	while ((1 << stages) < n) ++stages;
	mStages = stages;

	// 计算反转表
	mReverseTable.reserve(mN);
	for (int i = 0; i < mN; i ++)
	{
		mReverseTable.push_back(ReverseBitsLoop(i, mStages));
	}

}

void FFT::buildTwiddle()
{
    // 构造旋转矩阵，根据旋转矩阵的周期性，这里只需要存储n个数值就可以了
    const double PI = std::acos(-1.0);
    
    for (int i = 0; i < mN / 2; i ++)
    {
        double theta = -2.0 * PI * i / mN;
        mTw[i] = cd(cos(theta), sin(theta));
    }
}

std::vector<cd> FFT::forward(const std::vector<cd>& x) const
{
    // 1、根据原始序列位串反转表的计算，并将原始序列位串反转
	std::vector<cd> dst(mN);
	for (int i = 0; i < mN; i++)
	{
		dst[i] = x[mReverseTable[i]];
	}

	cd* data = dst.data();

	// 2、蝶形运算
	for (int s = 1; s <= mStages; ++s)
	{
		uint32_t blockLength = 1 << s; // 每一级的蝶形块的长度，例如，
										// 第1级就是相邻的两个元素做蝶形运算，那么蝶形块的长度就是2，第二级就是4，第三级就是8

		uint32_t halfBlockLength = blockLength >> 1;  // 每一级的蝶形块的长度的一半，也是两个蝶形运算输入数据的索引“距离”

		// 每一个块做蝶形运算
		for (uint32_t start = 0; start < mN; start += blockLength)
		{
			cd* low = data + start;
			cd* high = low + halfBlockLength;

			for (uint32_t j = 0; j < halfBlockLength; j ++)
			{
				const cd u = low[j];
				const cd v = high[j] * mTw[j * pow(2, mStages - s)];

				low[j] = u + v;
				high[j] = u - v;
			}
		}
	}
    return dst;
}

std::vector<cd> FFT::inverse(const std::vector<cd>& X) const
{
	// 1、根据原始序列位串反转表的计算，并将原始序列位串反转
	std::vector<cd> dst(mN);
	for (int i = 0; i < mN; i++)
	{
		dst[i] = X[mReverseTable[i]];
	}

	cd* data = dst.data();

	// 2、蝶形运算
	for (int s = 1; s <= mStages; ++s)
	{
		uint32_t blockLength = 1 << s; // 每一级的蝶形块的长度，例如，
		// 第1级就是相邻的两个元素做蝶形运算，那么蝶形块的长度就是2，第二级就是4，第三级就是8

		uint32_t halfBlockLength = blockLength >> 1;  // 每一级的蝶形块的长度的一半，也是两个蝶形运算输入数据的索引“距离”

		// 每一个块做蝶形运算
		for (uint32_t start = 0; start < mN; start += blockLength)
		{
			cd* low = data + start;
			cd* high = low + halfBlockLength;

			for (uint32_t j = 0; j < halfBlockLength; j++)
			{
				const cd u = low[j];
				const cd v = high[j] * std::conj(mTw[j * pow(2, mStages - s)]);

				low[j] = u + v;
				high[j] = u - v;
			}
		}
	}

	double divFactor = 1.0 / mN;
	for (int i = 0; i < mN; i++)
	{
		dst[i] *= divFactor;
	}

	return dst;
}

void TestFFT()
{
	std::vector<cd> input;
	input.push_back(cd(2.3));
	input.push_back(cd(3.0));
	input.push_back(cd(3.8));
	input.push_back(cd(4.2));
	input.push_back(cd(5.1));
	input.push_back(cd(4.5));
	input.push_back(cd(3.2));
	input.push_back(cd(1.8));

	FFT fft(input.size());
	std::vector<cd> output = fft.forward(input);

	printf("fft forward : \n");
	for (const auto& iter : output)
	{
		printf("real = %lf, imag = %lf\n", iter.real(), iter.imag());
	}

	std::vector<cd> inverseOutput = fft.inverse(output);

	printf("fft inverse : \n");
	for (const auto& iter : inverseOutput)
	{
		printf("real = %lf, imag = %lf\n", iter.real(), iter.imag());
	}
}
