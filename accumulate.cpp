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

#include "accumulate.h"

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
    phred_offset(decode_value_by_key< uint8_t >("phred offset", ontology)),
    platform(decode_value_by_key< Platform >("platform", ontology)),
    resolution(decode_value_by_key< int32_t >("resolution", ontology)),
    capacity(0),
    shortest(numeric_limits< int32_t >::max()),
    nucleic_acid_count_by_code(IUPAC_CODE_SIZE, 0) {

    } catch(ConfigurationError& error) {
        throw ConfigurationError("SegmentAccumulator :: " + error.message);

    } catch(exception& error) {
        throw InternalError("SegmentAccumulator :: " + string(error.what()));
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
        Value cycle_quality_report(kObjectType);
        Value cycle_nucleotide_quality_reports(kArrayType);
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
                    cycle_nucleotide_quality_reports.PushBack(cycle_nucleotide_quality_report, allocator);
                } else {
                    cycle_quality_report.AddMember("cycle quality distribution", cycle_quality_distribution, allocator);
                }
            }
        }
        container.AddMember("cycle nucleotide quality reports", cycle_nucleotide_quality_reports, allocator);
        container.AddMember("cycle nucleotide quality report", cycle_quality_report, allocator);

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

/*  ChannelAccumulator */

ChannelAccumulator::ChannelAccumulator(const Value& ontology) try :
    index(decode_value_by_key< int32_t >("index", ontology)),
    algorithm(decode_value_by_key< Algorithm >("algorithm", ontology)),
    disable_quality_control(decode_value_by_key< bool >("disable quality control", ontology)),
    undetermined(decode_value_by_key< bool >("undetermined", ontology)),
    concentration(decode_value_by_key< double >("concentration", ontology)),
    barcode(ontology),
    rg(ontology),
    count(0),
    multiplex_distance(0),
    multiplex_confidence(0),
    pf_count(0),
    pf_multiplex_distance(0),
    pf_multiplex_confidence(0),
    pf_fraction(0),
    pooled_fraction(0),
    pf_pooled_fraction(0),
    pooled_multiplex_fraction(0),
    pf_pooled_multiplex_fraction(0),
    accumulated_multiplex_distance(0),
    accumulated_multiplex_confidence(0),
    accumulated_pf_multiplex_distance(0),
    accumulated_pf_multiplex_confidence(0),
    segment_by_index(decode_value_by_key< vector< SegmentAccumulator > >("output feed by segment", ontology)) {

    } catch(ConfigurationError& error) {
        throw ConfigurationError("ChannelAccumulator :: " + error.message);

    } catch(exception& error) {
        throw InternalError("ChannelAccumulator :: " + string(error.what()));
};
void ChannelAccumulator::finalize(const OutputAccumulator& decoder_accumulator) {
    if(count > 0) {
        multiplex_distance = accumulated_multiplex_distance / double(count);
        multiplex_confidence = accumulated_multiplex_confidence / double(count);
        pf_fraction = double(pf_count) / double(count);
    }
    if(pf_count > 0) {
        pf_multiplex_distance = accumulated_pf_multiplex_distance / double(pf_count);
        pf_multiplex_confidence = accumulated_pf_multiplex_confidence / double(pf_count);
    }
    if(decoder_accumulator.count > 0) {
        pooled_fraction = double(count) / double(decoder_accumulator.count);
    }
    if(decoder_accumulator.pf_count > 0) {
        pf_pooled_fraction = double(pf_count) / double(decoder_accumulator.pf_count);
    }
    if(decoder_accumulator.multiplex_count > 0) {
        pooled_multiplex_fraction = double(count) / double(decoder_accumulator.multiplex_count);
    }
    if(decoder_accumulator.pf_multiplex_count > 0) {
        pf_pooled_multiplex_fraction = double(pf_count) / double(decoder_accumulator.pf_multiplex_count);
    }
    for(auto& accumulator : segment_by_index) {
        accumulator.finalize();
    }
};
ChannelAccumulator& ChannelAccumulator::operator+=(const ChannelAccumulator& rhs) {
    count += rhs.count;
    pf_count += rhs.pf_count;
    accumulated_multiplex_distance += rhs.accumulated_multiplex_distance;
    accumulated_multiplex_confidence += rhs.accumulated_multiplex_confidence;
    accumulated_pf_multiplex_distance += rhs.accumulated_pf_multiplex_distance;
    accumulated_pf_multiplex_confidence += rhs.accumulated_pf_multiplex_confidence;
    for(size_t i(0); i < segment_by_index.size(); ++i) {
        segment_by_index[i] += rhs.segment_by_index[i];
    }
    return *this;
};
template<> vector< ChannelAccumulator > decode_value_by_key(const Value::Ch* key, const Value& container) {
    vector< ChannelAccumulator > value;
    Value::ConstMemberIterator reference = container.FindMember(key);
    if(reference != container.MemberEnd()) {
        if(!reference->value.IsNull()) {
            if(reference->value.IsObject()) {
                value.reserve(reference->value.MemberCount());
                for(auto& record : reference->value.GetObject()) {
                    value.emplace_back(record.value);
                }
            } else { throw ConfigurationError(string(key) + " element must be a dictionary"); }
        }
    }
    return value;
};
bool encode_value(const ChannelAccumulator& value, Value& container, Document& document) {
    if(container.IsObject()) {
        encode_key_value("index", value.index, container, document);
        encode_value(value.rg, container, document);
        if(value.is_not_undetermined()) {
            encode_key_value("concentration", value.concentration, container, document);
            encode_key_value("barcode", value.barcode, container, document);
        }
        encode_key_value("count", value.count, container, document);
        if(value.is_not_undetermined()) {
            encode_key_value("multiplex distance", value.multiplex_distance, container, document);
            if(value.algorithm == Algorithm::PAMLD) {
                encode_key_value("multiplex confidence", value.multiplex_confidence, container, document);
            }
        }
        encode_key_value("pf count", value.pf_count, container, document);
        if(value.is_not_undetermined()) {
            encode_key_value("pf multiplex distance", value.pf_multiplex_distance, container, document);
            if(value.algorithm == Algorithm::PAMLD) {
                encode_key_value("pf multiplex confidence", value.pf_multiplex_confidence, container, document);
            }
        }
        encode_key_value("pf fraction", value.pf_fraction, container, document);
        encode_key_value("pooled fraction", value.pooled_fraction, container, document);
        encode_key_value("pf pooled fraction", value.pf_pooled_fraction, container, document);
        if(value.is_not_undetermined()) {
            encode_key_value("pooled multiplex fraction", value.pooled_multiplex_fraction, container, document);
            encode_key_value("pf pooled multiplex fraction", value.pf_pooled_multiplex_fraction, container, document);
        }
        if(!value.disable_quality_control) {
            Value feed_report_array(kArrayType);
            for(auto& accumulator : value.segment_by_index) {
                Value feed_report(kObjectType);
                encode_value(accumulator, feed_report, document);
                feed_report_array.PushBack(feed_report.Move(), document.GetAllocator());
            }
            container.AddMember("segment quality reports", feed_report_array.Move(), document.GetAllocator());
        }
        return true;
    } else { throw ConfigurationError("channel accumulator element must be a dictionary"); }
};
bool encode_key_value(const string& key, const ChannelAccumulator& value, Value& container, Document& document) {
    if(container.IsObject()) {
        Value element(kObjectType);
        encode_value(value, element, document);
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), element, document.GetAllocator());
        return true;
    } else { throw ConfigurationError(key + " element must be a dictionary"); }
};

/*  InputAccumulator */

InputAccumulator::InputAccumulator(const Value& ontology) try :
    disable_quality_control(decode_value_by_key< bool >("disable quality control", ontology)),
    count(0),
    pf_count(0),
    pf_fraction(0),
    segment_by_index(decode_value_by_key< vector< SegmentAccumulator > >("input feed by segment", ontology)) {

    } catch(ConfigurationError& error) {
        throw ConfigurationError("InputAccumulator :: " + error.message);

    } catch(exception& error) {
        throw InternalError("InputAccumulator :: " + string(error.what()));
};
void InputAccumulator::finalize() {
    if(count > 0) {
        pf_fraction = double(pf_count) / double(count);
    }
    for(auto& accumulator : segment_by_index) {
        accumulator.finalize();
    }
};
InputAccumulator& InputAccumulator::operator+=(const InputAccumulator& rhs) {
    count += rhs.count;
    pf_count += rhs.pf_count;

    for(size_t i(0); i < segment_by_index.size(); ++i) {
        segment_by_index[i] += rhs.segment_by_index[i];
    }
    return *this;
};
bool encode_key_value(const string& key, const InputAccumulator& value, Value& container, Document& document) {
    if(container.IsObject()) {
        Value element(kObjectType);
        encode_key_value("count", value.count, element, document);
        encode_key_value("pf count", value.pf_count, element, document);
        encode_key_value("pf fraction", value.pf_fraction, element, document);
        if(!value.disable_quality_control) {
            Value feed_report_array(kArrayType);
            for(auto& accumulator : value.segment_by_index) {
                Value feed_report(kObjectType);
                encode_value(accumulator, feed_report, document);
                feed_report_array.PushBack(feed_report.Move(), document.GetAllocator());
            }
            element.AddMember("segment quality reports", feed_report_array.Move(), document.GetAllocator());
        }
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), element.Move(), document.GetAllocator());
        return true;
    } else { throw ConfigurationError(key + " element must be a dictionary"); }
};

/*  OutputAccumulator */

OutputAccumulator::OutputAccumulator(const Value& ontology) try :
    algorithm(decode_value_by_key< Algorithm >("algorithm", ontology)),
    count(0),
    multiplex_count(0),
    multiplex_fraction(0),
    multiplex_distance(0),
    multiplex_confidence(0),
    pf_count(0),
    pf_fraction(0),
    pf_multiplex_count(0),
    pf_multiplex_fraction(0),
    pf_multiplex_distance(0),
    pf_multiplex_confidence(0),
    multiplex_pf_fraction(0),
    accumulated_multiplex_distance(0),
    accumulated_multiplex_confidence(0),
    accumulated_pf_multiplex_distance(0),
    accumulated_pf_multiplex_confidence(0),
    channel_by_index(decode_value_by_key< vector< ChannelAccumulator > >("codec", ontology)),
    undetermined(find_value_by_key("undetermined", ontology)) {

    } catch(ConfigurationError& error) {
        throw ConfigurationError("OutputAccumulator :: " + error.message);

    } catch(exception& error) {
        throw InternalError("OutputAccumulator :: " + string(error.what()));
};
void OutputAccumulator::finalize() {
    count += undetermined.count;
    pf_count += undetermined.pf_count;

    for(auto& accumulator : channel_by_index) {
        count += accumulator.count;
        multiplex_count += accumulator.count;
        accumulated_multiplex_distance += accumulator.accumulated_multiplex_distance;
        if(algorithm == Algorithm::PAMLD) {
            accumulated_multiplex_confidence += accumulator.accumulated_multiplex_confidence;
        }
        pf_count += accumulator.pf_count;
        pf_multiplex_count += accumulator.pf_count;
        accumulated_pf_multiplex_distance += accumulator.accumulated_pf_multiplex_distance;
        if(algorithm == Algorithm::PAMLD) {
            accumulated_pf_multiplex_confidence += accumulator.accumulated_pf_multiplex_confidence;
        }
    }
    if(count > 0) {
        multiplex_fraction = double(multiplex_count) / double(count);
        multiplex_distance = accumulated_multiplex_distance / double(count);
        if(algorithm == Algorithm::PAMLD) {
            multiplex_confidence = accumulated_multiplex_confidence / double(count);
        }
        pf_fraction = double(pf_count) / double(count);
    }
    if(pf_count > 0) {
        pf_multiplex_fraction = double(pf_multiplex_count) / double(pf_count);
        pf_multiplex_distance = accumulated_pf_multiplex_distance / double(pf_count);
        if(algorithm == Algorithm::PAMLD) {
            pf_multiplex_confidence = accumulated_pf_multiplex_confidence / double(pf_count);
        }
    }
    if(multiplex_count > 0) {
        multiplex_pf_fraction = double(pf_multiplex_count) / double(multiplex_count);
    }
    undetermined.finalize(*this);
    for(auto& accumulator : channel_by_index) {
        accumulator.finalize(*this);
    }

};
OutputAccumulator& OutputAccumulator::operator+=(const OutputAccumulator& rhs) {
    undetermined += rhs.undetermined;
    for(size_t i(0); i < channel_by_index.size(); ++i) {
        channel_by_index[i] += rhs.channel_by_index[i];
    }
    return *this;
};
bool encode_key_value(const string& key, const OutputAccumulator& value, Value& container, Document& document) {
    if(container.IsObject()) {
        Value element(kObjectType);
        encode_key_value("count", value.count, element, document);
        encode_key_value("multiplex count", value.multiplex_count, element, document);
        encode_key_value("multiplex fraction", value.multiplex_fraction, element, document);
        encode_key_value("multiplex distance", value.multiplex_distance, element, document);
        if(value.algorithm == Algorithm::PAMLD) {
            encode_key_value("multiplex confidence", value.multiplex_confidence, element, document);
        }
        encode_key_value("pf count", value.pf_count, element, document);
        encode_key_value("pf fraction", value.pf_fraction, element, document);
        encode_key_value("pf multiplex count", value.pf_multiplex_count, element, document);
        encode_key_value("pf multiplex fraction", value.pf_multiplex_fraction, element, document);
        encode_key_value("pf multiplex distance", value.pf_multiplex_distance, element, document);
        if(value.algorithm == Algorithm::PAMLD) {
            encode_key_value("pf multiplex confidence", value.pf_multiplex_confidence, element, document);
        }
        encode_key_value("multiplex pf fraction", value.multiplex_pf_fraction, element, document);
        encode_key_value("undetermined quality report", value.undetermined, element, document);

        Value channel_report_array(kArrayType);
        for(auto& channel : value.channel_by_index) {
            Value channel_report(kObjectType);
            encode_value(channel, channel_report, document);
            channel_report_array.PushBack(channel_report.Move(), document.GetAllocator());
        }
        element.AddMember("read group quality reports", channel_report_array.Move(), document.GetAllocator());

        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), element.Move(), document.GetAllocator());
        return true;
    } else { throw ConfigurationError(key + " container is not a dictionary"); }
};
