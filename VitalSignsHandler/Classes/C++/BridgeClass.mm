#import "BridgeClass.h"
#import "Patch.hpp"
@interface BridgeClass ()
@property (nonatomic)Mircod_Patch_class* helper;
@end

@implementation BridgeClass

-(instancetype)init{
    if (self = [super init]){
        self.helper = new Mircod_Patch_class;
    }
    return self;
}

-(NSDictionary*)updateECG:(int)data{
    ECG_data ecg =  self.helper->ECG_signal_class.update_ECG(data);
    NSMutableDictionary *dic = [NSMutableDictionary new];
    [dic setValue:[NSNumber numberWithDouble:ecg.ECG] forKey:@"ECG"];
    [dic setValue:[NSNumber numberWithDouble:ecg.AvgHR] forKey:@"ECG_HR"];
    [dic setValue:[NSNumber numberWithDouble:ecg.AvgRR] forKey:@"ECG_RRInterval"];
    [dic setValue:[NSNumber numberWithDouble:ecg.HRV] forKey:@"ECG_HRV"];
    [dic setValue:[NSNumber numberWithDouble:ecg.Resp1] forKey:@"ECG_BreathingRate"];
    return [dic copy];
    
}


-(NSDictionary*)updatePPGRed:(double)Red IR:(double)IR{
    SPO2_data ppg = self.helper->PPG_signal_class.update_PPG(Red,IR);
    NSMutableDictionary *dic = [NSMutableDictionary new];
    [dic setValue:[NSNumber numberWithDouble:ppg.SpO2] forKey:@"SPO2"];
    [dic setValue:[NSNumber numberWithDouble:ppg.AvgHR] forKey:@"SPO2_HR"];
    [dic setValue:[NSNumber numberWithDouble:ppg.HRV] forKey:@"SPO2_HRV"];
    return [dic copy];
}

@end
