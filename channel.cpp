/*
    Pheniqs : PHilology ENcoder wIth Quality Statistics
    Copyright (C) 2018  Lior Galanti
    NYU Center for Genetics and System Biology

    Author: Lior Galanti <lior.galanti@nyu.edu>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "channel.h"

/*  NucleotideAccumulator */

NucleotideAccumulator::NucleotideAccumulator() :
    count(0),
    min_quality(0),
    max_quality(0),
    sum_quality(0),
    mean_quality(0),
    Q1(0),
    Q3(0),
    IQR(0),
    LW(0),
    RW(0),
    median_quality(0),
    distribution(EFFECTIVE_PHRED_RANGE, 0) {
};
void NucleotideAccumulator::finalize() {
    for(auto& q : distribution) {
        count += q;
    }
    if(count > 0) {
        for(size_t q(0); q < distribution.size(); ++q) {
            const uint64_t value(distribution[q]);
            sum_quality += (value * q);
            if(value != 0) {
                max_quality = q;
                if(min_quality == 0) {
                    min_quality = q;
                }
            }
        }
        mean_quality = double(sum_quality) / double(count);
        median_quality = quantile(0.5);
        Q1 = quantile(0.25);
        Q3 = quantile(0.75);
        IQR = Q3 - Q1;

        double W(Q1 - IQR * 1.5);
        LW = (W < min_quality) ? min_quality : W;

        W = Q3 + IQR * 1.5;
        RW = (W > max_quality) ? max_quality : W;
    }
};
NucleotideAccumulator& NucleotideAccumulator::operator=(const NucleotideAccumulator& rhs) {
    if(this != &rhs) {
        count = rhs.count;
        min_quality = rhs.min_quality;
        max_quality = rhs.max_quality;
        sum_quality = rhs.sum_quality;
        mean_quality = rhs.mean_quality;
        Q1 = rhs.Q1;
        Q3 = rhs.Q3;
        IQR = rhs.IQR;
        LW = rhs.LW;
        RW = rhs.RW;
        median_quality = rhs.median_quality;
        distribution = rhs.distribution;
    }
    return *this;
};
NucleotideAccumulator& NucleotideAccumulator::operator+=(const NucleotideAccumulator& rhs) {
    for(size_t q(0); q < distribution.size(); ++q) {
        distribution[q] += rhs.distribution[q];
    }
    return *this;
};

/*  CycleAccumulator */

CycleAccumulator::CycleAccumulator() :
    nucleotide_by_code(IUPAC_CODE_SIZE) {
};
void CycleAccumulator::finalize() {
    /* accumulate all nucleotide variations in the NO_NUCLEOTIDE accumulative distribution */
    for(uint8_t i(1); i < nucleotide_by_code.size(); ++i) {
        for(uint8_t p(0); p < EFFECTIVE_PHRED_RANGE; ++p) {
            nucleotide_by_code[NO_NUCLEOTIDE].distribution[p] += nucleotide_by_code[i].distribution[p];
        }
    }
    for(auto& distribution : nucleotide_by_code) {
        distribution.finalize();
    }
};
CycleAccumulator& CycleAccumulator::operator=(const CycleAccumulator& rhs) {
    if(this != &rhs) {
        nucleotide_by_code = rhs.nucleotide_by_code;
    }
    return *this;
};
CycleAccumulator& CycleAccumulator::operator+=(const CycleAccumulator& rhs) {
    for(size_t i(0); i < nucleotide_by_code.size(); ++i) {
        nucleotide_by_code[i] += rhs.nucleotide_by_code[i];
    }
    return *this;
};

/*  AveragePhreadAccumulator */

AveragePhreadAccumulator::AveragePhreadAccumulator() :
    count(0),
    min_value(0),
    max_value(0),
    sum_value(0),
    mean_value(0),
    distribution(EFFECTIVE_PHRED_RANGE, 0) {
};
void AveragePhreadAccumulator::finalize() {
    if(count > 0) {
        mean_value = sum_value / double(count);
    }
};
AveragePhreadAccumulator& AveragePhreadAccumulator::operator=(const AveragePhreadAccumulator& rhs) {
    if(this != &rhs) {
        count = rhs.count;
        min_value = rhs.min_value;
        max_value = rhs.max_value;
        sum_value = rhs.sum_value;
        mean_value = rhs.mean_value;
        distribution = rhs.distribution;
    }
    return *this;
};
AveragePhreadAccumulator& AveragePhreadAccumulator::operator+=(const AveragePhreadAccumulator& rhs) {
    count += rhs.count;
    sum_value += rhs.sum_value;
    min_value = min(min_value, rhs.min_value);
    max_value = max(max_value, rhs.max_value);
    for(size_t i(0); i < distribution.size(); ++i) {
        distribution[i] += rhs.distribution[i];
    }
    return *this;
};

/*  SegmentAccumulator */

SegmentAccumulator::SegmentAccumulator(const Value& ontology) try :
    index(decode_value_by_key< int32_t >("index", ontology)),
    url(decode_value_by_key< URL >("url", ontology)),
    capacity(0),
    shortest(numeric_limits< int32_t >::max()),
    nucleic_acid_count_by_code(IUPAC_CODE_SIZE, 0) {

    } catch(Error& error) {
        error.push("SegmentAccumulator");
        throw;
};
void SegmentAccumulator::finalize() {
    if(shortest == numeric_limits< int32_t >::max()) {
        shortest = 0;
    }
    for(auto& c : cycle_by_index) {
        c.finalize();
    }
    average_phred.finalize();
};
SegmentAccumulator& SegmentAccumulator::operator+=(const SegmentAccumulator& rhs) {
    if(rhs.capacity > capacity) {
        cycle_by_index.resize(rhs.capacity);
        capacity = rhs.capacity;
    }
    shortest = min(shortest, rhs.shortest);
    for(uint8_t c(0); c < nucleic_acid_count_by_code.size(); ++c) {
        nucleic_acid_count_by_code[c] += rhs.nucleic_acid_count_by_code[c];
    }

    for(int32_t i(0); i < rhs.capacity; ++i) {
        cycle_by_index[i] += rhs.cycle_by_index[i];
    }
    average_phred += rhs.average_phred;
    return *this;
};
template<> vector< SegmentAccumulator > decode_value_by_key(const Value::Ch* key, const Value& container) {
    if(!container.IsNull()) {
        Value::ConstMemberIterator reference = container.FindMember(key);
        if(reference != container.MemberEnd()) {
            if(!reference->value.IsNull()) {
                vector< SegmentAccumulator > value;
                if(reference->value.IsArray() && !reference->value.Empty()) {
                    value.reserve(reference->value.Size());
                    for(const auto& element : reference->value.GetArray()) {
                        value.emplace_back(element);
                    }
                }
                return value;
            } else { throw ConfigurationError(string(key) + " is null"); }
        } else { throw ConfigurationError(string(key) + " not found"); }
    } else { throw ConfigurationError(string(key) + " container is null"); }
};
bool encode_value(const SegmentAccumulator& value, Value& container, Document& document) {
    if(container.IsObject()) {
        Document::AllocatorType& allocator = document.GetAllocator();

        encode_key_value("url", value.url, container, document);
        encode_key_value("min sequence length", value.shortest, container, document);
        encode_key_value("max sequence length", value.capacity, container, document);
        Value quality_control_by_cycle(kObjectType);
        Value quality_control_by_nucleotide(kArrayType);
        for(uint8_t n(0); n < value.nucleic_acid_count_by_code.size(); ++n) {
            if(value.nucleic_acid_count_by_code[n] > 0) {
                Value cycle_quality_distribution(kObjectType);
                Value cycle_count(kArrayType);
                Value cycle_quality_first_quartile(kArrayType);
                Value cycle_quality_third_quartile(kArrayType);
                Value cycle_quality_interquartile_range(kArrayType);
                Value cycle_quality_left_whisker(kArrayType);
                Value cycle_quality_right_whisker(kArrayType);
                Value cycle_quality_min(kArrayType);
                Value cycle_quality_max(kArrayType);
                Value cycle_quality_mean(kArrayType);
                Value cycle_quality_median(kArrayType);
                for(size_t c(0); c < value.cycle_by_index.size(); ++c) {
                    cycle_count.PushBack(Value(value.cycle_by_index[c].nucleotide_by_code[n].count).Move(), allocator);
                    cycle_quality_first_quartile.PushBack(Value(value.cycle_by_index[c].nucleotide_by_code[n].Q1).Move(), allocator);
                    cycle_quality_third_quartile.PushBack(Value(value.cycle_by_index[c].nucleotide_by_code[n].Q3).Move(), allocator);
                    cycle_quality_interquartile_range.PushBack(Value(value.cycle_by_index[c].nucleotide_by_code[n].IQR).Move(), allocator);
                    cycle_quality_left_whisker.PushBack(Value(value.cycle_by_index[c].nucleotide_by_code[n].LW).Move(), allocator);
                    cycle_quality_right_whisker.PushBack(Value(value.cycle_by_index[c].nucleotide_by_code[n].RW).Move(), allocator);
                    cycle_quality_min.PushBack(Value(value.cycle_by_index[c].nucleotide_by_code[n].min_quality).Move(), allocator);
                    cycle_quality_max.PushBack(Value(value.cycle_by_index[c].nucleotide_by_code[n].max_quality).Move(), allocator);
                    cycle_quality_mean.PushBack(Value(value.cycle_by_index[c].nucleotide_by_code[n].mean_quality).Move(), allocator);
                    cycle_quality_median.PushBack(Value(value.cycle_by_index[c].nucleotide_by_code[n].median_quality).Move(), allocator);
                }
                cycle_quality_distribution.AddMember("cycle count", cycle_count, allocator);
                cycle_quality_distribution.AddMember("cycle quality first quartile", cycle_quality_first_quartile, allocator);
                cycle_quality_distribution.AddMember("cycle quality third quartile", cycle_quality_third_quartile, allocator);
                cycle_quality_distribution.AddMember("cycle quality interquartile range", cycle_quality_interquartile_range, allocator);
                cycle_quality_distribution.AddMember("cycle quality left whisker", cycle_quality_left_whisker, allocator);
                cycle_quality_distribution.AddMember("cycle quality right whisker", cycle_quality_right_whisker, allocator);
                cycle_quality_distribution.AddMember("cycle quality min", cycle_quality_min, allocator);
                cycle_quality_distribution.AddMember("cycle quality max", cycle_quality_max, allocator);
                cycle_quality_distribution.AddMember("cycle quality mean", cycle_quality_mean, allocator);
                cycle_quality_distribution.AddMember("cycle quality median", cycle_quality_median, allocator);
                if(n > 0) {
                    Value cycle_nucleotide_quality_report(kObjectType);
                    encode_key_value("nucleotide count", value.nucleic_acid_count_by_code[n], cycle_nucleotide_quality_report, document);
                    encode_key_value("nucleotide", string(1, BamToAmbiguousAscii[n]), cycle_nucleotide_quality_report, document);
                    cycle_nucleotide_quality_report.AddMember("cycle quality distribution", cycle_quality_distribution, allocator);
                    quality_control_by_nucleotide.PushBack(cycle_nucleotide_quality_report, allocator);
                } else {
                    quality_control_by_cycle.AddMember("cycle quality distribution", cycle_quality_distribution, allocator);
                }
            }
        }
        container.AddMember("quality control by nucleotide", quality_control_by_nucleotide, allocator);
        container.AddMember("quality control by cycle", quality_control_by_cycle, allocator);

        Value average_phred_report(kObjectType);
        encode_key_value("average phred score min", value.average_phred.min_value, average_phred_report, document);
        encode_key_value("average phred score max", value.average_phred.max_value, average_phred_report, document);
        encode_key_value("average phred score mean", value.average_phred.mean_value, average_phred_report, document);

        Value segment_accumulator(kArrayType);
        for(size_t i(0); i < value.average_phred.distribution.size(); ++i) {
            segment_accumulator.PushBack(Value(value.average_phred.distribution[i]).Move(), allocator);
        }
        average_phred_report.AddMember("average phred score distribution", segment_accumulator, allocator);
        container.AddMember("average phred score report", average_phred_report, allocator);
        return true;
    } else { throw ConfigurationError("feed accumulator element must be a dictionary"); }
};

/*  Channel */

Channel::Channel(const Value& ontology) try :
    Barcode(ontology),
    rg(ontology),
    include_filtered(decode_value_by_key< bool >("include filtered", ontology)),
    enable_quality_control(decode_value_by_key< bool >("enable quality control", ontology)),
    output_feed_url_by_segment(decode_value_by_key< list< URL > >("output", ontology)),
    segment_accumulator_by_index(decode_value_by_key< vector< SegmentAccumulator > >("output feed by segment", ontology["feed"])) {

    } catch(Error& error) {
        error.push("Channel");
        throw;
};
Channel::Channel(const Channel& other) :
    Barcode(other),
    rg(other.rg),
    include_filtered(other.include_filtered),
    enable_quality_control(other.enable_quality_control),
    output_feed_url_by_segment(other.output_feed_url_by_segment),
    output_feed_lock_order(other.output_feed_lock_order),
    output_feed_by_segment(other.output_feed_by_segment) {
};
void Channel::populate(unordered_map< URL, Feed* >& feed_by_url) {
    map< int32_t, Feed* > feed_by_index;

    /* populate the output feed by segment array */
    output_feed_by_segment.reserve(output_feed_url_by_segment.size());
    for(const auto& url : output_feed_url_by_segment) {
        Feed* feed(feed_by_url[url]);
        output_feed_by_segment.emplace_back(feed);
        if(feed_by_index.count(feed->index) == 0) {
            feed_by_index.emplace(make_pair(feed->index, feed));
        }
    }
    output_feed_by_segment.shrink_to_fit();

    /* populate the output feed lock order array */
    output_feed_lock_order.reserve(feed_by_index.size());
    for(auto& record : feed_by_index) {
        /* /dev/null is not really being written to so we don't need to lock it */
        if(!record.second->is_dev_null()) {
            output_feed_lock_order.push_back(record.second);
        }
    }
    output_feed_lock_order.shrink_to_fit();
};
void Channel::finalize(const AccumulatingDecoder& parent) {
    AccumulatingIdentifier::finalize(parent);
    if(enable_quality_control) {
        for(auto& segment_accumulator : segment_accumulator_by_index) {
            segment_accumulator.finalize();
        }
    }
};
void Channel::encode(Value& container, Document& document) const {
    Barcode::encode(container, document);
    if(container.IsObject()) {
        encode_value(rg, container, document);
        if(enable_quality_control) {
            Value quality_control_by_segment(kArrayType);
            for(auto& accumulator : segment_accumulator_by_index) {
                Value feed_report(kObjectType);
                encode_value(accumulator, feed_report, document);
                quality_control_by_segment.PushBack(feed_report.Move(), document.GetAllocator());
            }
            container.AddMember("quality control by segment", quality_control_by_segment.Move(), document.GetAllocator());
        }
    } else { throw ConfigurationError("element must be a dictionary"); }
};
Channel& Channel::operator+=(const Channel& rhs) {
    Barcode::operator+=(rhs);
    if(enable_quality_control) {
        for(size_t index(0); index < segment_accumulator_by_index.size(); ++index) {
            segment_accumulator_by_index[index] += rhs.segment_accumulator_by_index[index];
        }
    }
    return *this;
};
template<> vector< Channel > decode_value_by_key(const Value::Ch* key, const Value& container) {
    vector< Channel > value;
    Value::ConstMemberIterator reference = container.FindMember(key);
    if(reference != container.MemberEnd()) {
        value.reserve(reference->value.MemberCount());
        for(auto& record : reference->value.GetObject()) {
            value.emplace_back(record.value);
        }
    }
    return value;
};
